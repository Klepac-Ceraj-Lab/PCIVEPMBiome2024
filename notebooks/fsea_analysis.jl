#####
# FSEA
#####

function perform_fsea(infile::String, outfile::String; should_consolidate=false)

    df = subset(CSV.read(infile, DataFrame), "pvalue" => ByRow(!isnan))
    Ts = df.z
    # clamp extreme |z| to avoid undue leverage
    Ts = map(x -> isfinite(x) ? clamp(x, -7, 7) : NaN, Ts)

    genes = map(f -> replace(f, "UniRef90_" => ""), df.feature)
    neuroactive = Leap.getneuroactive(genes; consolidate=should_consolidate)


    results = ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        if isempty(ixs) || length(ixs) < 10
            return (; geneset=gs, enrichment=NaN, pvalue=NaN)
        end
        valid = findall(!isnan, Ts)
        ixs = intersect(ixs, valid)
        if isempty(ixs)
            return (; geneset=gs, enrichment=NaN, pvalue=NaN)
        end

        fes = fsea(MWU(), Ts, ixs)
        (; geneset=gs,
           enrichment = enrichment_score(fes.setranks, fes.nfeatures),
           pvalue = PCIVEPMBiome2024.pvalue(fes))
    end

    tmp = DataFrame(results)
    subset!(tmp, :pvalue => ByRow(!isnan))
    tmp.qvalue = adjust(tmp.pvalue, BenjaminiHochberg())
    sort!(tmp, [:qvalue, :pvalue])
    CSV.write(outfile, tmp)
    tmp
    
end

## Model 1
perform_fsea( "manuscript/FSEA/lms_Model1_mbiome_age.csv", "manuscript/FSEA/fsea_full_Model1_mbiome_age.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model1_mbiome_age.csv", "manuscript/FSEA/fsea_consolidated_Model1_mbiome_age.csv"; should_consolidate=true)

perform_fsea( "manuscript/FSEA/lms_Model1_child_sex.csv", "manuscript/FSEA/fsea_full_Model1_child_sex.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model1_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model1_child_sex.csv"; should_consolidate=true)

## Model 2
perform_fsea( "manuscript/FSEA/lms_Model2_child_sex.csv", "manuscript/FSEA/fsea_full_Model2_child_sex.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model2_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model2_child_sex.csv"; should_consolidate=true)

perform_fsea( "manuscript/FSEA/lms_Model2_visual.csv", "manuscript/FSEA/fsea_full_Model2_visual.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model2_visual.csv", "manuscript/FSEA/fsea_consolidated_Model2_visual.csv"; should_consolidate=true)

perform_fsea( "manuscript/FSEA/lms_Model2_entropy.csv", "manuscript/FSEA/fsea_full_Model2_entropy.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model2_entropy.csv", "manuscript/FSEA/fsea_consolidated_Model2_entropy.csv"; should_consolidate=true)

## Model 3
perform_fsea( "manuscript/FSEA/lms_Model3_child_sex.csv", "manuscript/FSEA/fsea_full_Model3_child_sex.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model3_child_sex.csv", "manuscript/FSEA/fsea_consolidated_Model3_child_sex.csv"; should_consolidate=true)

perform_fsea( "manuscript/FSEA/lms_Model3_visual.csv", "manuscript/FSEA/fsea_full_Model3_visual.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model3_visual.csv", "manuscript/FSEA/fsea_consolidated_Model3_visual.csv"; should_consolidate=true)

perform_fsea( "manuscript/FSEA/lms_Model3_entropy.csv", "manuscript/FSEA/fsea_full_Model3_entropy.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model3_entropy.csv", "manuscript/FSEA/fsea_consolidated_Model3_entropy.csv"; should_consolidate=true)

perform_fsea( "manuscript/FSEA/lms_Model3_product.csv", "manuscript/FSEA/fsea_full_Model3_product.csv"; should_consolidate=false)
perform_fsea( "manuscript/FSEA/lms_Model3_product.csv", "manuscript/FSEA/fsea_consolidated_Model3_product.csv"; should_consolidate=true)