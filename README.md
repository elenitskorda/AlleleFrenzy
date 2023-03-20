
# Author: Eleni Theofania Skorda
# Date:   21/03/2023
# Course: BINP29 DNA Sequencing Analysis

# Title: **Allele Frenzy** - A simulator of change in minor allele frequencies over time in multiple populations with different initial allele frequencies and SNP values


## Description

>This is a Shiny app that generates changes of minor allele frequencies of population over different generations. The app allows users to define different parameters such as number of populations, number of generations, number of SNPs, the starting allele frequency, and the rate of drift. The app then creates a plot of the minor allele frequency for each SNP over time.

## Packages

>For the execution of the app, you will need the shiny and the ggplot libraries witht the below code. All the dependencies are installed with them.

```r
install.packages("shiny")
install.packages("ggplot2")
```
> **shiny**  In order to generate the interface for a web application, we will use the shiny packages which is amongst the most popular popular packages for this purpose. Thus it will make it easier for the user to generate **the minor allele frequencies (MAF)** by using the second library which is the ggplot2

> **ggplot2** is a package for generating graphics/plots in R. It is built on a layered grammar of graphics and provides a flexible and powerful system for visualizing data. In this code, ggplot2 is used to create a plot of MAF changes over time (generations).

> Once you have the packages installed, you can run the app by running the shinyApp() function with the ui and server functions as arguments:

```r
library(shiny)
library(ggplot2)

shinyApp(ui, server)
```
## Usage 

> Once the app is running, you can adjust specific input parameters using **ONLY** the sliders. These parameters are:
1. Number of populations
2. Number of generations
3. Number of SNPs
4. Sstarting allele frequency
5. Rate of drift

 While changing the above values, the plot in the main panel called Variation Voyage will update to show the MAF for each SNP for each generation.

## Code Structure

> The code is divided into two parts: the user interface (UI) and the server. The UI defines the layout of the app and the input controls, while the server defines the logic for the simulation and the plot. The UI is defined using the fluidPage() function from the shiny package, which allows us to create a responsive, fluid layout that adjusts to the size of the user's browser window. The UI consists of a sidebar layout with input controls for the simulation parameters and a main panel with the output plot. The server is defined using the server() function, which contains the logic for the simulation and the plot. This is described by this part of the code:

```r

# Define user interface (UI) elements
ui <- fluidPage(
  # Define the overall page layout with a sidebar and main panel
  sidebarLayout(
     # Sidebar panel containing input controls for the simulation parameters
    sidebarPanel(
      # Title and parameter selection options for the simulation
      # Headers and numeric inputs
      h2("AlleleFrenzy"),
      # Headers and numeric inputs
      h3("Select Parameters"),
      # Numeric input for number of populations
      numericInput("num_pops", "Number of populations:", value = 4, min = 1),
      # Numeric input for number of generations
      numericInput("num_generations", "Number of generations:", value = 20, min = 1),
      # Numeric input for number of SNPs
      numericInput("num_snps", "Number of SNPs:", value = 5, min = 1),
      # Numeric input for Starting allele frequency
      numericInput("starting_af", "Starting allele frequency:", value = 0.5, min = 0, max = 1, step = 0.01),
      # Numeric input for genetic drift
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
```


>The server code is included in the observe() function, which allows us to update the plot in real time as the user changes the values of the parameters. The simulation itself is performed using a for loop that iterates over the number of generations specified by the user. In each iteration of the loop, the server introduces genetic drift to the allele frequencies, creates a new mixed population, calculates the minor allele frequency for each SNP, and saves the results to a data frame. 

```r
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

```
>Once the simulation is complete, the data frames are combined and plotted using ggplot2.  The output$minor_af_plot object is generated using renderPlot() function. The renderPlot() function takes a plot generation function (in this case, ggplot()) as its argument and generates a plot based on the output of that function. Specifically, the renderPlot() function is used to generate a plot of minor allele frequency changes over time, using the ggplot2 package which creates the plot. The plot is created based on a data frame called ''af_df'' .This data frame contains information on theMAF for each SNP over generations. This can be described by this part:

```r
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
```
>The plot need to be like the pictuce below:

![image](https://user-images.githubusercontent.com/32453911/226433913-d9cca898-dca2-4694-9920-f58d47c92f99.png)
>Figure 1: Demonstration of Allele Frenzy application tool. On the left side of the window, the user adjusts the input parameters, and the plot are adjust-ed accordingly to the right side with the generation change on the y axis and the minor allele frequency on the x axis.

## Limitations and Future Work

>There are some limitations to the current version of the Allele Frenzy app which are documented below:
1. Does not represent of real-world genetic traits.
2. Single gene locus at a time
3. Assume idealized conditions
>In order to improve and make the app more realistic, we need to consider polygenic inheritance and multiple gene loci simultaneously. Additionally, including more realistic population dynamics and evolutionary processes could make the app even more useful for examining genetic variation and evolution. Overall, the "Allele Frenzy" app is a useful tool for educational purposes, but there is still room for further development and improvement.

## Conclusion

>In conclusion, the "Allele Frenzy" Shiny app is a simple yet useful tool for visualizing the changes in MAFs over different generations in a mixed population of multiple populations. The app is very easy to use and gives a visual representation of the simulated genetic drift in the population. The ggplot2 package is used to create the plot, and the Shiny package allows for the user interface to be built and for the app to be deployed on the web. This app can be useful for researchers or students interested in population genetics, evolutionary biology, or genetic drift, and can also serve as an educational tool to teach these concepts to a broader audience.



