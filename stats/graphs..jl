using CSV
using DataFrames
using Plots

df = DataFrame(CSV.File("data.csv", delim = '\t'))
print(df)
