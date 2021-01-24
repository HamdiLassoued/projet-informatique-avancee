#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#Importation de librairies
library("shiny")
library("shinydashboard")
library("ggplot2")
library("deSolve")
library("cowplot")
library("ggrepel")
library("tidyverse")

#Le modele SIR
sir <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
        dS <- -beta * S * I
        dI <-  beta * S * I - gamma * I
        dR <-                 gamma * I
        dV <- 0
        return(list(c(dS, dI, dR, dV)))
    })
}

# UI

ui <- dashboardPage(
    dashboardHeader(title ='SIR avec vaccin'),
    dashboardSidebar(
        sliderInput(
            "connum",
            "Numero reproductif de base R0 ",
            min = .5,
            max = 20,
            value = 5
        ),
        sliderInput(
            "pinf",
            "Nombre infectes au debut",
            min = 1,
            max = 100,
            value = 10
        ),
        sliderInput(
            "pvac",
            "Proportion de vaccines + precedemment immunises (%)",
            min = 0,
            max = 100,
            value = 75
        ),
        sliderInput(
            "vaceff",
            "Efficacite du vaccin (%)",
            min = 0,
            max = 100,
            value = 85
        ),
        sliderInput(
            "infper",
            "Periode infection (jours)",
            min = 1,
            max = 30,
            value = 7
        ),
        sliderInput(
            "timeframe",
            "Periode (jours)",
            min = 1,
            max = 500,
            value = 200
        )    
    ),
    dashboardBody(
        tags$head(tags$style(HTML('
                              /* main sidebar 
                              .skin-blue .main-sidebar {
                              background-color: #808080;
                              } */

                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #fffff8;
                              }                              
                              '))),
        
        fluidRow(plotOutput("distPlot")),
        br(),
        fluidRow(
            # definir des valueBoxes dynamiques
            valueBoxOutput("progressBox", width = 6),
            valueBoxOutput("approvalBox", width = 6),
            valueBoxOutput("BRRBox", width = 6),
            valueBoxOutput("HIBox", width = 6)
            
        ),
        br(),
        br()
    )
)

# preparation pour le plot
server <- function(input, output) {
    # input reactive
    popsize <- 11690000   # population tunisienne
    dataInput <- reactive({
        init       <-
            c(
                S = 1 - input$pinf / popsize - input$pvac / 100 * input$vaceff / 100,
                I = input$pinf / popsize,
                R = 0,
                V = input$pvac / 100 * input$vaceff / 100
            )
        ## beta: parametre de transmission et gamma: parametre de guerison
        parameters <-
            c(beta = input$connum * 1 / input$infper,
              (1 - input$pvac/100*input$vaceff/100),
              gamma = 1 / input$infper)
        ## temps frame
        times <- seq(0, input$timeframe, by = .2)
        
        ## Resolution avec ode
        out <-
            ode(
                y = init,
                times = times,
                func = sir,
                parms = parameters
            )
        as.data.frame(out)
    })
    
    output$distPlot <- renderPlot({
        out <-
            dataInput() %>%
            gather(key, value, -time) %>%
            mutate(
                id = row_number(),
                key2 = recode(
                    key,
                    S = "Susceptibles (S)",
                    I = "Infectes (I)",
                    R = "Retablies (R)",
                    V = "Vaccines + precedemment immunises"
                ),
                keyleft = recode(
                    key,
                    S = "Susceptibles (S)",
                    I = "",
                    R = "",
                    V = "Vaccines + precedemment immunises"
                ),
                keyright = recode(
                    key,
                    S = "",
                    I = "Infectes (I)",
                    R = "Retablies (R)",
                    V = ""
                )
            )
        
        ggplot(data = out,
               aes(
                   x = time,
                   y = value,
                   group = key2,
                   col = key2,
                   label = key2,
                   data_id = id
               )) +  ylim(0, 1) +
            ylab("Pourcentage de la population totale") + xlab("Temps (jours)") +
            geom_line(size = 2) +
            geom_text_repel(
                data = subset(out, time == max(time)),
                aes(label = keyright),
                size = 6,
                segment.size  = 0.2,
                segment.color = "grey",
                nudge_x = 0,
                hjust = 1,
                direction = "y"
            ) +
            geom_text_repel(
                data = subset(out, time == min(time)),
                aes(label = keyleft),
                size = 6,
                segment.size  = 0.2,
                segment.color = "grey",
                nudge_x = 0,
                hjust = 0,
                direction = "y"
            ) +
            theme(legend.position = "none") +
            scale_colour_manual(values = c("red", "green", "black", "blue")) +
            scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
            theme(
                rect=element_rect(size=0),
                legend.position="none",
                panel.background=element_rect(fill="transparent", colour=NA),
                plot.background =element_rect(fill="transparent", colour=NA),
                legend.key =     element_rect(fill = "transparent", colour = "transparent")
            )
        
    })
    
    output$progressBox <- renderValueBox({
        valueBox(
            dataInput() %>% filter(time == max(time)) %>% select(R) %>% mutate(R = round(100 *
                                                                                             R, 2)) %>% paste0("%"),
            "Proportion de la population totale touchee par la maladie",
            icon = icon("thermometer"),
            color = "red"
        )
    })
    
    output$approvalBox <- renderValueBox({
        valueBox(
            paste0(round(100 * (dataInput() %>% filter(row_number() == n()) %>% mutate(res = (R + I) / (S + I + R)) %>% pull("res") ), 2), "%"),
            "Proportion de personnes sensibles a la maladie",
            icon = icon("exclamation"),
            color = "orange"
        )
    })
    
    output$BRRBox <- renderValueBox({
        valueBox(
            paste0(round(input$connum *(1 - input$pvac / 100 * input$vaceff / 100), 2), ""),
            "R0 effective (pour la population au debut de l'epidemie, lorsque l'immunite est prise en compte)",
            icon = icon("arrow-circle-right"),
            color = "green"
        )
    })
    
    # Temps manquants
    output$HIBox <- renderValueBox({
        valueBox(
            paste0(round(100 * (1 - 1 / (input$connum) ), 2), "%"),
            "Proportion de la population qui doit etre immunisee pour obtenir l immunite du troupeau au debut d une epidemie",
            icon = icon("stethoscope"),
            color = "blue"
        )
    })
    
}

# Execution
shinyApp(ui = ui, server = server)
