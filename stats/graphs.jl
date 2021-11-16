using CSV
using DataFrames
using Plots

df = DataFrame(CSV.File("data.csv", delim = '\t'));

# Divide all the values by 1000 to turn them into ms.
for i = 2:6
  df[!, i] = df[!, i] ./ 1000;
end

bar([1 2 3], df.comYoutube)

