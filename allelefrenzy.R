# Author: Eleni Theofania Skorda
# Author: Eleni Theofania Skorda
# Date:   21/03/2023
# Course: BINP29 DNA Sequencing Analysis

# Title: **Allele Frenzy** - A simulator of change in minor allele frequencies over time in multiple populations with different initial allele frequencies and SNP values


#The script creates a Shiny web application that simulates the changes in minor allele frequencies
#of a set of single nucleotide polymorphisms (SNPs) across multiple populations over a specified 
#number of generations. The application allows users to adjust the simulation parameters, including
#the number of populations, number of generations, number of SNPs, starting allele frequency, and 
#the degree of genetic drift.

#The user interface (UI) consists of a sidebar and a main panel. The sidebar contains numeric input
#controls for the simulation parameters, and the main panel displays the simulation output plot. 
#The UI is created using the fluidPage() function from the shiny package.

#The server logic is implemented in the server() function. The function uses the observe() function
#to listen for changes in the input parameters and generate new output accordingly.The server
#function starts by generating the initial allele frequencies for each population using the matrix()
#function. It then initializes an empty list to store the allele frequency data over time. Next,
#it enters a for loop that simulates the changes in allele frequencies over the specified number 
#of generations.

#In each iteration of the loop, the function introduces genetic drift using the rnorm() function,
#which generates random numbers from a normal distribution. It then ensures that the allele 
#frequencies remain between 0 and 1. The function then creates a mixed population by taking a 
#weighted average of the allele frequencies across all populations for each SNP.

#The function then calculates the minor allele frequency (MAF) for each SNP using the apply() 
#function, which applies a function to each column of a matrix. It saves the MAF data for the 
#current generation to a data frame and appends it to the list of allele frequencies over time.

#After the loop completes, the function combines the data frames into one using the do.call() 
#and rbind() functions. Finally, it creates a plot of the MAF changes over time using the ggplot2 
#package. The plot shows a line graph of the MAF for each SNP over the specified number of
#generations, with each SNP in a separate subplot.


# Load necessary libraries
library(shiny)
library(ggplot2)

# Define user interface (UI) elements
ui <- fluidPage(
  # Define the overall page layout with a sidebar and main panel
  sidebarLayout(
    # Sidebar panel containing input controls for the simulation parameters
    sidebarPanel(
      # Title and parameter selection options for the simulation
      # Headers and numeric inputs
      h2("AlleleFrenzy"),
      h3("Select Parameters"),
      # Numeric input for number of populations
      numericInput("num_pops", "Number of populations:", value = 4, min = 1),
      # Numeric input for number of generations
      numericInput("num_generations", "Number of generations:", value = 20, min = 1),
      # Numeric input for number of SNPs
      numericInput("num_snps", "Number of SNPs:", value = 5, min = 1),
      # Numeric input for Starting allele frequency
      numericInput("starting_af", "Starting allele frequency:", value = 0.5, min = 0, max = 1, step = 0.01),
      # Numeric input for genetic drift rate
      numericInput("drift", "Drift:", value = 0.01, min = 0, max = 1, step = 0.01)
      
    ),
    # Main panel containing the simulation output plot
    mainPanel(
      # Title and plot for the simulation output
      h3("Variation Voyage"),
      plotOutput("minor_af_plot")
    )
  )
)

server <- function(input, output) {
  observe({
    # Generate the initial allele frequencies
    af <- matrix(input$starting_af, nrow = input$num_pops, ncol = input$num_snps)
    
    # Initialize an empty list to store allele frequencies over time
    af_list <- list()
    
    # Simulate allele frequency changes over time
    for (i in 1:input$num_generations) {
      # Introduce genetic drift
      af <- af + rnorm(input$num_pops*input$num_snps, 0, input$drift)
      # Ensure allele frequencies remain between 0 and 1
      af[af < 0] <- 0
      af[af > 1] <- 1
      
      # Create new mixed population
      mixed_af <- matrix(NA, nrow = 1, ncol = input$num_snps)
      for (j in 1:input$num_snps) {
        # Calculate weighted average of allele frequencies
        mixed_af[1, j] <- weighted.mean(af[, j], sample(input$num_pops, input$num_pops, replace = TRUE))
      }
      
      
      # Calculate minor allele frequency for all SNPs
      minor_af <- apply(af, 2, function(x) min(x + 0.01, 1-x + 0.01))
      
      # Save minor allele frequencies to data frame
      af_df <- data.frame(
        generation = i,
        population = "Mixed",
        snp = 1:input$num_snps,
        af = minor_af  # only the minor allele frequency
      )
      
      # Append data frame to list
      af_list[[i]] <- af_df
      
    }
    
    # Combine data frames into one data frame
    af_df <- do.call(rbind, af_list)
    
    # Plot minor allele frequency changes over time
    output$minor_af_plot <- renderPlot({
      # Create a ggplot with af_df as the data frame and define the aesthetics
      ggplot(af_df, aes(x = generation, y = af, color = factor(snp))) +
        # Add a line plot for each SNP showing its minor allele frequency over time
        geom_line() +
        # Add axis and legend labels
        labs(x = "Generation", y = "Minor Allele Frequency", color = "SNP") +
        # Add facetting to separate the plots by SNP
        facet_wrap(~snp, ncol = 3) +
        # Keep the color scale consistent across all plots
        scale_color_discrete(drop = FALSE)
    })
  })
}

shinyApp(ui, server)