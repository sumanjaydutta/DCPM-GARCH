# Install necessary packages if not installed
if (!require("httr")) install.packages("httr")
if (!require("rvest")) install.packages("rvest")
if (!require("quantmod")) install.packages("quantmod")

# Load libraries
library(httr)
library(rvest)
library(quantmod)

# Step 1: Scrape Nifty 500 Stock Names using httr
nifty_500_url <- "https://www1.nseindia.com/live_market/dynaContent/live_watch/stock_watch/nifty500StockWatch.htm"

# Use httr to get page content
response <- GET(nifty_500_url, user_agent("Mozilla/5.0"))

# Check if the connection was successful
if (status_code(response) == 200) {
  
  # Parse the HTML content
  nifty_500_page <- read_html(content(response, "text"))
  
  # Extract the stock names (Modify this based on structure of the table)
  nifty_500_table <- nifty_500_page %>%
    html_node("table") %>%
    html_table(fill = TRUE)
  
  # Extract first column (stock symbols or names)
  stock_names <- nifty_500_table[[1]]  # Modify as per the table's structure
  stock_names <- stock_names %>% as.character()
  
  print(stock_names)  # Check the list of stock names
  
} else {
  stop("Failed to connect to the NSE website.")
}

# Step 2: Fetch Daily Closing Prices using quantmod
# Define the start and end dates
start_date <- as.Date("2022-01-01")
end_date <- as.Date("2023-01-01")

# Create a list to store the closing prices
closing_prices_list <- list()

# Loop through each stock and fetch the daily closing prices
for (stock in stock_names) {
  stock_symbol <- paste0(stock, ".NS")  # NSE stock symbols in Yahoo
  
  # Try fetching the stock data, adding a delay to avoid rate-limiting issues
  tryCatch({
    Sys.sleep(2)  # Add a 2-second delay between requests
    stock_data <- getSymbols(stock_symbol, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
    
    # Extract closing prices
    closing_prices <- Cl(stock_data)
    
    # Store the closing prices in the list
    closing_prices_list[[stock]] <- closing_prices
  }, error = function(e) {
    message(paste("Could not retrieve data for stock:", stock_symbol))
  })
}

# View the fetched closing prices
print(closing_prices_list)

# Save to CSV if needed
# saveRDS(closing_prices_list, "nifty500_closing_prices.rds")
