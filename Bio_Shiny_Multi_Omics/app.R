# shiny_app/app.R

# Load required libraries
library(shiny)
library(ggplot2)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Data Analysis - Bioinformatics Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      selectInput("x_var", "Select X Variable", choices = NULL),
      selectInput("y_var", "Select Y Variable", choices = NULL)
    ),
    mainPanel(
      plotOutput("scatterPlot"),
      tableOutput("summaryTable")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to load data
  data <- reactive({
    req(input$file)  # Ensure file is uploaded
    read.csv(input$file$datapath)
  })
  
  # Update UI elements based on the uploaded data
  observe({
    df <- data()
    updateSelectInput(session, "x_var", choices = colnames(df))
    updateSelectInput(session, "y_var", choices = colnames(df))
  })
  
  # Render the scatter plot based on selected variables
  output$scatterPlot <- renderPlot({
    df <- data()
    ggplot(df, aes_string(x = input$x_var, y = input$y_var)) +
      geom_point() +
      theme_minimal()
  })
  
  # Render the summary table based on the uploaded data
  output$summaryTable <- renderTable({
    df <- data()
    summary(df)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
