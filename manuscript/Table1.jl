mdata_df = select(khula_pci_mbiome_data, mdata_cols)

# Initialize an empty DataFrame to store the combined summary data
combined_summary = DataFrame(
    :column => "Data points",
    :value => "Count",
    :count => nrow(mdata_df),
    :percentage => "100%"
)

# Summarize each column and calculate percentages
sex_summary = DataFrames.combine(groupby(mdata_df, :child_sex), nrow => :count)
sex_summary.percentage = 100 * sex_summary[:, :count] ./ nrow(mdata_df)
sex_summary.column .= "Child sex"
rename!(sex_summary, :child_sex => :value)
sort!(sex_summary, :value)
    
# Append the summary to the combined DataFrame
combined_summary = vcat(combined_summary, sex_summary)

## Dealing with EPDS
age_data = mdata_df[:, "mbiome_sample_age"]
age_min_max = DataFrame(
    column = "Age",
    value = "Median [Min, Max]",
    count = median(age_data),
    percentage = "[$(minimum(age_data)), $(maximum(age_data))]"
)

mated_sd = mdata_df[:, "MaternalEducation"]
mated_mean_sd = DataFrame(
    column = "Maternal Education (formal years)",
    value = "Mean (SD)",
    count = mean(mated_sd),
    percentage = std(mated_sd)
)

matent_sd = mdata_df[:, "MaternalEntropy"]
matent_mean_sd = DataFrame(
    column = "Maternal Entropy",
    value = "Mean (SD)",
    count = mean(matent_sd),
    percentage = std(matent_sd)
)

infvis_sd = mdata_df[:, "InfantVisAtt"]
infvis_mean_sd = DataFrame(
    column = "Visual Attention Shifts",
    value = "Mean (SD)",
    count = mean(infvis_sd),
    percentage = std(infvis_sd)
)

combined_summary = vcat(combined_summary, age_min_max, mated_mean_sd, matent_mean_sd, infvis_mean_sd)

# Display the combined summary table
using PrettyTables
pretty_table(combined_summary)
CSV.write("manuscript/tables/Table1.csv", combined_summary)