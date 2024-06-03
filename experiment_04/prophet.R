
library(prophet)

df <- read.csv('https://raw.githubusercontent.com/facebook/prophet/main/examples/example_wp_log_peyton_manning.csv')

m <- prophet(df)

future <- make_future_dataframe(m, periods = 365)
tail(future)

forecast <- predict(m, future)
tail(forecast[c('ds', 'yhat', 'yhat_lower', 'yhat_upper')])

plot(m, forecast)

prophet_plot_components(m, forecast)


# trend flexibility 
# default = 0.05
# higher = more flexible
m <- prophet(df, changepoint.prior.scale = 0.5)
forecast <- predict(m, future)
plot(m, forecast)


# fourier order of seasonalities
# default = 10
# higher = more flexible (fits higher-freq changes)
m <- prophet(df, yearly.seasonality = 20)
prophet:::plot_yearly(m)


# different seasonality freq
# ???
m <- prophet(weekly.seasonality=FALSE)
m <- add_seasonality(m, name='monthly', period=30.5, fourier.order=5)
m <- fit.prophet(m, df)
forecast <- predict(m, future)
prophet_plot_components(m, forecast)


# additional regressors SEE SECTION
nfl_sunday <- function(ds) {
  dates <- as.Date(ds)
  month <- as.numeric(format(dates, '%m'))
  as.numeric((weekdays(dates) == "Sunday") & (month > 8 | month < 2))
}
df$nfl_sunday <- nfl_sunday(df$ds)

m <- prophet()
m <- add_regressor(m, 'nfl_sunday')
m <- fit.prophet(m, df)

future$nfl_sunday <- nfl_sunday(future$ds)

forecast <- predict(m, future)
plot(m, forecast)
prophet_plot_components(m, forecast)


# monthly data
df <- read.csv('https://raw.githubusercontent.com/facebook/prophet/main/examples/example_retail_sales.csv')
m <- prophet(df, seasonality.mode = 'multiplicative')
future <- make_future_dataframe(m, periods = 120, freq = 'month')
fcst <- predict(m, future)
plot(m, fcst)
prophet_plot_components(m, fcst)

# In monthly data, yearly seasonality can also be modeled with binary extra 
# regressors. In particular, the model can use 12 extra regressors like 
# is_jan, is_feb, etc. where is_jan is 1 if the date is in Jan and 0 otherwise. 
# This approach would avoid the within-month unidentifiability seen above. 
# Be sure to use yearly_seasonality=False if monthly extra regressors are being added.
