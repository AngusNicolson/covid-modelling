
# Load Packages
library(shiny)
library(deSolve)

ui <- fluidPage(
    # Set mathjax running
    withMathJax(),
    # Application title
    titlePanel("COVID-19 Compartmental Models"),
 
    sidebarLayout(
        sidebarPanel(
            h2("Simulation Parameters"),
            fluidRow(column(6,
                            sliderInput("time.max",
                                        "Number of days:",
                                        min = 1,
                                        max = 365,
                                        value = 62)),
                     column(6,
                            numericInput("initial.S",
                                         "Susceptible",
                                         min=1e2,
                                         max=1e8,
                                         value=10000))),
            fluidRow(column(4,
                            numericInput("initial.Ii",
                                         "Incubation",
                                         min=0,
                                         max=1e8,
                                         value=10)),
                     column(4,
                            numericInput("initial.Ip",
                                         "Pre-symptomatic",
                                         min=0,
                                         max=1e8,
                                         value=10)),
                     column(4,
                            numericInput("initial.Is",
                                         "Symptomatic",
                                         min=0,
                                         max=1e8,
                                         value=10))),
            h2("Model Options"),
            checkboxInput("lockdown", "Lockdown"),
            conditionalPanel(
                condition = "input.lockdown == true",
                sliderInput("lockdown.range", "Lockdown Days", min = 0, 
                            max = 365, value = c(30, 50)),
                sliderInput("lockdown.value", "Lockdown Effectiveness",
                             min=0, max=1, value=0.9)
            ),
            checkboxInput("hospital", "Hospital Capacity"),
            conditionalPanel(
                condition = "input.hospital == true",
                numericInput("hospital.capacity", "Max Hospital Capacity",
                            min=0, max=1e8, value=10000),
                numericInput("hospital.value", "Over-Capacity Penalty",
                            min=0, max=100, value=5)
            ),
            h2("Model Parameters"),
            fluidRow(column(4,
                            numericInput("beta_p",
                                         "$$\\beta_p$$",
                                         min=0,
                                         max=10,
                                         value=0.2)),
                     column(4,
                            numericInput("beta_s",
                                         "$$\\beta_s$$",
                                         min=0,
                                         max=10,
                                         value=0.55))),
            fluidRow(column(4,
                            numericInput("k_ip",
                                         "$$k_{ip}$$",
                                         min=0,
                                         max=10,
                                         value=0.5)),
                     column(4,
                            numericInput("k_ps",
                                         "$$k_{ps}$$",
                                         min=0,
                                         max=10,
                                         value=0.5)),
                     column(4,
                            numericInput("k_sr",
                                         "$$k_{sr}$$",
                                         min=0,
                                         max=10,
                                         value=0.2))),
            fluidRow(column(4,
                            numericInput("w",
                                         "$$w$$",
                                         min=0,
                                         max=10,
                                         value=0.002)),
                     column(4,
                            numericInput("delta_Is",
                                         "$$\\delta_{I_s}$$",
                                         min=0,
                                         max=10,
                                         value=0.001))),
            width=5
        ),

        mainPanel(
            tabsetPanel(
                type = "tabs",
                tabPanel("Plot", 
                         fluidRow(
                             column(6, plotOutput("allPlot")),
                             column(6, plotOutput("totalPlot"))
                         ),
                         fluidRow(
                             column(6, plotOutput("rPlot")),
                             column(6, plotOutput("incidencePlot"))
                         )
                ),
                tabPanel("Diagram", 
                         imageOutput("diagram")
                )
            ),
            width=7
        )
    )
)

server <- function(input, output) {
    
    output$infectedPlot <- renderPlot({
        out <- get_model_output(input)
        
        plot(out$time, out$Is, type = "l", col = "red", lwd = 2,
             xlab = "Time", ylab = "Number", xaxs="i", yaxs="i")
        lines(out$time, out$Ii, lwd = 2, col = "yellow")
        lines(out$time, out$Ip, lwd = 2, col = "orange")
        legend("topleft", legend = c("Incubation", "Pre-Symptomatic", "Symptomatic"), lty = 1, col = c("orange", "yellow", "red"), lwd = 2, bty = "n")
    })
    
    output$incidencePlot <- renderPlot({
        out <- get_model_output(input)
        N   <- rowSums(out[,2:6])
        
        incidence <- (input$beta_p * out$Ip + input$beta_s * out$Is) * out$S / N
        death_incidence <- input$delta_Is * out$Is
        
        #plot(out$time, incidence, type = "l", col = "red", lwd = 2,
        #     xlab = "Time", ylab = "No. Cases a day", xaxs="i", yaxs="i")
        #lines(out$time, death_incidence, lwd = 2, col = "black")
        #legend("topleft", legend = c("Infected", "Deaths"), lty = 1, col = c("red", "black"), lwd = 2, bty = "n")
        
        par(mar=c(5, 4, 4, 6) + 0.1)
        
        ## Plot first set of data and draw its axis
        plot(out$time, incidence, type = "l", col = "red", lwd = 2, axes=FALSE, ylim=c(0,max(incidence)), xlab="", ylab="")
        axis(2, ylim=c(0,max(incidence)),col="black")  ## las=1 makes horizontal labels
        mtext("No. new infections a day",side=2,line=2.5)
        box()
        
        ## Allow a second plot on the same graph
        par(new=TRUE)
        
        ## Plot the second plot and put axis scale on right
        plot(out$time, death_incidence, type="l",  xlab="", ylab="", ylim=c(0,max(death_incidence)), axes=FALSE, col="black", lwd=2)
        ## a little farther out (line=4) to make room for labels
        mtext("No. deaths a day",side=4,col="black",line=4) 
        axis(4, ylim=c(0,max(death_incidence)), col="black",col.axis="black")
        
        ## Draw the time axis
        axis(1,pretty(range(out$time),10))
        mtext("Time / Days",side=1,col="black",line=2.5)  
        
        ## Add Legend
        legend("topleft",legend=c("Infections","Deaths"), lty=1, lwd=2, col=c("red", "black"), bty="n")
        
    })
    
    output$totalPlot <- renderPlot({
        out <- get_model_output(input)
        I <- rowSums(out[,3:5])
        
        par(mar=c(5, 4, 4, 6) + 0.1)
        
        ## Plot first set of data and draw its axis
        plot(out$time, I, type = "l", col = "red", lwd = 2, axes=FALSE, ylim=c(0,max(I)), xlab="", ylab="")
        axis(2, ylim=c(0,max(I)),col="black")  ## las=1 makes horizontal labels
        mtext("No. infected",side=2,line=2.5)
        lines(out$time, out$Ii, lty=1, lwd = 1.5, col = "#ffb366")  # 30, 60, 100 
        lines(out$time, out$Ip, lty=1, lwd = 1.5, col = "#ff8040") # 20, 75, 100
        lines(out$time, out$Is, lty=1, lwd = 1.5, col = "#ff3f19") # 10, 90, 100
        
        box()
        
        ## Allow a second plot on the same graph
        par(new=TRUE)
        
        ## Plot the second plot and put axis scale on right
        plot(out$time, out$D, type="l",  xlab="", ylab="", ylim=c(0,max(out$D)), axes=FALSE, col="black", lwd=2)
        ## a little farther out (line=4) to make room for labels
        mtext("No. deaths",side=4,col="black",line=4) 
        axis(4, ylim=c(0,max(out$D)), col="black",col.axis="black")
        
        ## Draw the time axis
        axis(1,pretty(range(out$time),10))
        mtext("Time / Days",side=1,col="black",line=2.5)  
        
        ## Add Legend
        legend("topleft",legend=c("Infections","Deaths"), lty=1, lwd=2, col=c("red", "black"), bty="n")
        legend("left", legend = c("Incubation", "Pre-Symptomatic", "Symptomatic"), lty = 1, col = c("#ffb366", "#ff8040", "#ff3f19"), lwd = 1, bty = "n")
    })
    
    output$allPlot <- renderPlot({
        out <- get_model_output(input)
        N   <- rowSums(out[,2:6])
        I   <- rowSums(out[,3:5])
        plot(out$time, out$S, type = "l", col = "blue", lwd = 2, ylim = c(0, (max(N))),
             xlab = "Time", ylab = "Number", xaxs="i", yaxs="i")
        lines(out$time, I, lwd = 2, col = "red")
        lines(out$time, out$R, lwd = 2, col = "green")
        lines(out$time, out$D, lwd = 2, col = "black")
        
        legend("left", legend = c("Susceptible", "Infected", 'Recovered', "Dead"), lty = 1, col = c("blue", "red", "green", "black"), lwd = 2, bty = "n")
    })
    
    output$rPlot <- renderPlot({
        out <- get_model_output(input)
        N   <- rowSums(out[,2:6])
        r0  <- input$beta_p / input$k_ps + input$beta_s / (input$delta_Is + input$k_sr)
        r   <- (out[, "S"] / N) * r0
        
        if (input$lockdown) {
            lockdown.min <- input$lockdown.range[[1]]
            lockdown.max <- input$lockdown.range[[2]]
            
            lockdown.mask <- seq(0, input$time.max, by = 1)
            lockdown.mask <- (lockdown.mask < lockdown.max) & (lockdown.mask > lockdown.min)
            normal.mask   <- 1 - lockdown.mask
            lockdown.mask <- lockdown.mask * (1 - input$lockdown.value)
            mask          <- normal.mask + lockdown.mask
            r <- r * mask
        }
        
        plot(out[, "time"], r, type="l", xaxs="i", yaxs="i", ylim=c(0, max(r)), xlab = "Time", ylab = "R")
        abline(h=1, lty=2)
        text(max(out[, "time"])*0.90, r0*0.95, paste("R0: ", format(r0, digits=3)))
    })
    
    output$deathsPlot <- renderPlot({
        out <- get_model_output(input)
        plot(out[, "time"], out[, "D"], type="l", xaxs="i", yaxs="i" , xlab = "Time", ylab = "Deaths")
    })
    
    output$diagram <- renderImage({
        return(list(
            src = "images/COVID-19 SIRS model.png",
            contentType = "image/png",
            alt = "COVID-19 SIRS Model Diagram"
        ))
    }, deleteFile = FALSE)
}

get_model_output <- function(input) {
    
    if (input$lockdown) {
        lockdown.min <- input$lockdown.range[[1]]
        lockdown.max <- input$lockdown.range[[2]]
    } else {
        lockdown.min <- 0
        lockdown.max <- 0
    }
    
    ps <- c(beta_p = input$beta_p, beta_s = input$beta_s,
            k_ip = input$k_ip, k_ps = input$k_ps, k_sr=input$k_sr,
            delta_Is = input$delta_Is, w=input$w,
            lockdown.min=lockdown.min, lockdown.max=lockdown.max, lockdown.value=input$lockdown.value
            )
    state <- c(S = input$initial.S, Ii = input$initial.Ii, Ip=input$initial.Ip, Is = input$initial.Is, R=0, D=0)
    times <- seq(0, input$time.max, by = 1)
    
    out <- ode(y = state, times = times, func = SIRS_model, parms = ps)
    out = as.data.frame(out)
    out
}

SIRS_model <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        if ((t > lockdown.min) && (t < lockdown.max) ){
            beta_p = beta_p * (1 - lockdown.value)
            beta_s = beta_s * (1 - lockdown.value)
        }
        N <- S + Ii + Ip + Is + R
        dS <- w * R - (beta_p * Ip + beta_s * Is) * S / N
        dIi <- (beta_p * Ip + beta_s * Is) * S / N - k_ip * Ii
        dIp <- k_ip * Ii - k_ps * Ip
        dIs <- k_ps * Ip - (k_sr + delta_Is) * Is
        dR <- k_sr * Is - w * R
        dD <- delta_Is * Is
        
        list(c(dS, dIi, dIp, dIs, dR, dD))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
