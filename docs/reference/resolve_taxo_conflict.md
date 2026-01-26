# Resolve taxonomic conflict in the tax_table of a phyloseq object

Resolve taxonomic conflict in the tax_table of a phyloseq object

## Usage

``` r
resolve_taxo_conflict(
  physeq,
  pattern_tax_ranks = NULL,
  method = c("consensus", "rel_majority", "abs_majority", "preference", "unanimity"),
  strict = FALSE,
  second_method = c("consensus", "rel_majority", "abs_majority", "preference",
    "unanimity"),
  nb_agree_threshold = 1,
  keep_tax_ranks = TRUE,
  new_names = NULL,
  preference_pattern = NULL,
  collapse_string = "/",
  replace_collapsed_rank_by_NA = FALSE
)
```

## Arguments

- physeq:

  (required) A phyloseq object.

- pattern_tax_ranks:

  (character vector, default NULL) A vector of patterns to aggregate
  taxonomic ranks. For example "^Genus" stands for all taxonomic ranks
  (columns in tax_table slot) starting with Genus.

- method:

  (character, default "consensus") One of "consensus", "rel_majority",
  "abs_majority", "preference" or "unanimity". See details in the
  documentation of the function
  [`MiscMetabar::resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.html).

- strict:

  (logical, default FALSE) If TRUE, NA are considered as informative in
  resolving conflict (i.e. NA are taken into account in vote). See
  details for more information.

- second_method:

  (character, default "consensus") One of "consensus", "rel_majority",
  "abs_majority", or "unanimity". Only used if method = "preference".
  See details.

- nb_agree_threshold:

  (integer, default 1) The minimum number of times a value must arise to
  be selected using vote. If 2, only taxonomic values present at least 2
  times in the vector are kept.

- keep_tax_ranks:

  (logical, default TRUE) If TRUE, keep the old taxonomic ranks in the
  result.

- new_names:

  (character vector, default NULL) A vector of new names for the
  taxonomic columns.

- preference_pattern:

  (character, default NULL) A pattern to match the only column used as
  preferred one if method = "preference".

- collapse_string:

  (character, default "/") The character to collapse taxonomic names
  when multiple assignment is done.

- replace_collapsed_rank_by_NA:

  (logical, default FALSE) If TRUE, all multiple assignments (all
  taxonomic rank including the 'collapse_string' parameter) are replaced
  by NA.

## Value

A phyloseq object

## Author

Adrien Taudière

## Examples

``` r
data_fungi_mini_new <- assign_sintax(data_fungi_mini,
  ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  ),
  behavior = "add_to_phyloseq"
)

resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "consensus", new_names = "Genus_consensus")@tax_table
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           Genus_consensus                               
#> ASV7   NA                  NA                                            
#> ASV8   NA                  "Stereum"                                     
#> ASV12  NA                  "Xylodon"                                     
#> ASV18  NA                  "Stereum"                                     
#> ASV25  NA                  "Ossicaulis/Tricholoma"                       
#> ASV26  NA                  "Stereum"                                     
#> ASV27  NA                  "Antrodiella"                                 
#> ASV29  NA                  "Basidiodendron"                              
#> ASV32  "Sistotrema_sp"     "Sistotrema"                                  
#> ASV34  NA                  "Entocybe"                                    
#> ASV35  NA                  "Fomes"                                       
#> ASV41  "Mycena_sp"         "Mycena"                                      
#> ASV42  NA                  "Ossicaulis"                                  
#> ASV46  NA                  "Helicogloea"                                 
#> ASV47  NA                  "Radulomyces"                                 
#> ASV48  NA                  "Elmerina"                                    
#> ASV49  NA                  "Phanerochaete"                               
#> ASV50  NA                  "Gloeohypochnicium"                           
#> ASV53  NA                  "Fomes"                                       
#> ASV54  NA                  "Auricularia"                                 
#> ASV58  NA                  "Fomes"                                       
#> ASV59  NA                  "Hyphoderma"                                  
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                                  
#> ASV62  "Serendipita_sp"    "Serendipita"                                 
#> ASV63  NA                  NA                                            
#> ASV64  NA                  "Trametes"                                    
#> ASV67  NA                  "Xylodon"                                     
#> ASV68  "Tricholoma_sp"     "Ossicaulis/Tricholoma"                       
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"              
#> ASV72  NA                  NA                                            
#> ASV75  "Peniophora_reidii" "Peniophora"                                  
#> ASV77  NA                  "Ossicaulis/Tricholoma"                       
#> ASV82  NA                  "Exidia"                                      
#> ASV83  NA                  NA                                            
#> ASV85  "Polyporales_sp"    "Peniophorella/Polyporales_gen_Incertae_sedis"
#> ASV91  NA                  "Auricularia"                                 
#> ASV93  NA                  NA                                            
#> ASV94  NA                  "Stereum"                                     
#> ASV99  NA                  "Fomes"                                       
#> ASV100 NA                  NA                                            
#> ASV101 NA                  "Marchandiomyces"                             
#> ASV104 NA                  "Hericium"                                    
#> ASV105 NA                  "Xylodon"                                     
#> ASV107 NA                  "Xylodon"                                     
#> ASV108 NA                  "Exidia"                                      
resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "consensus", new_names = "Genus_consensus", replace_collapsed_rank_by_NA = TRUE)@tax_table
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           Genus_consensus                 
#> ASV7   NA                  NA                              
#> ASV8   NA                  "Stereum"                       
#> ASV12  NA                  "Xylodon"                       
#> ASV18  NA                  "Stereum"                       
#> ASV25  NA                  NA                              
#> ASV26  NA                  "Stereum"                       
#> ASV27  NA                  "Antrodiella"                   
#> ASV29  NA                  "Basidiodendron"                
#> ASV32  "Sistotrema_sp"     "Sistotrema"                    
#> ASV34  NA                  "Entocybe"                      
#> ASV35  NA                  "Fomes"                         
#> ASV41  "Mycena_sp"         "Mycena"                        
#> ASV42  NA                  "Ossicaulis"                    
#> ASV46  NA                  "Helicogloea"                   
#> ASV47  NA                  "Radulomyces"                   
#> ASV48  NA                  "Elmerina"                      
#> ASV49  NA                  "Phanerochaete"                 
#> ASV50  NA                  "Gloeohypochnicium"             
#> ASV53  NA                  "Fomes"                         
#> ASV54  NA                  "Auricularia"                   
#> ASV58  NA                  "Fomes"                         
#> ASV59  NA                  "Hyphoderma"                    
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                    
#> ASV62  "Serendipita_sp"    "Serendipita"                   
#> ASV63  NA                  NA                              
#> ASV64  NA                  "Trametes"                      
#> ASV67  NA                  "Xylodon"                       
#> ASV68  "Tricholoma_sp"     NA                              
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                  NA                              
#> ASV75  "Peniophora_reidii" "Peniophora"                    
#> ASV77  NA                  NA                              
#> ASV82  NA                  "Exidia"                        
#> ASV83  NA                  NA                              
#> ASV85  "Polyporales_sp"    NA                              
#> ASV91  NA                  "Auricularia"                   
#> ASV93  NA                  NA                              
#> ASV94  NA                  "Stereum"                       
#> ASV99  NA                  "Fomes"                         
#> ASV100 NA                  NA                              
#> ASV101 NA                  "Marchandiomyces"               
#> ASV104 NA                  "Hericium"                      
#> ASV105 NA                  "Xylodon"                       
#> ASV107 NA                  "Xylodon"                       
#> ASV108 NA                  "Exidia"                        

resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "preference", preference_pattern = ".x$")@tax_table
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           ^Genus\\._preference            
#> ASV7   NA                  NA                              
#> ASV8   NA                  "Stereum"                       
#> ASV12  NA                  "Xylodon"                       
#> ASV18  NA                  "Stereum"                       
#> ASV25  NA                  "Ossicaulis"                    
#> ASV26  NA                  "Stereum"                       
#> ASV27  NA                  "Antrodiella"                   
#> ASV29  NA                  "Basidiodendron"                
#> ASV32  "Sistotrema_sp"     "Sistotrema"                    
#> ASV34  NA                  "Entocybe"                      
#> ASV35  NA                  "Fomes"                         
#> ASV41  "Mycena_sp"         "Mycena"                        
#> ASV42  NA                  "Ossicaulis"                    
#> ASV46  NA                  "Helicogloea"                   
#> ASV47  NA                  "Radulomyces"                   
#> ASV48  NA                  "Elmerina"                      
#> ASV49  NA                  "Phanerochaete"                 
#> ASV50  NA                  "Gloeohypochnicium"             
#> ASV53  NA                  "Fomes"                         
#> ASV54  NA                  "Auricularia"                   
#> ASV58  NA                  "Fomes"                         
#> ASV59  NA                  "Hyphoderma"                    
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                    
#> ASV62  "Serendipita_sp"    "Serendipita"                   
#> ASV63  NA                  NA                              
#> ASV64  NA                  "Trametes"                      
#> ASV67  NA                  "Xylodon"                       
#> ASV68  "Tricholoma_sp"     "Ossicaulis"                    
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                  NA                              
#> ASV75  "Peniophora_reidii" "Peniophora"                    
#> ASV77  NA                  "Ossicaulis"                    
#> ASV82  NA                  "Exidia"                        
#> ASV83  NA                  NA                              
#> ASV85  "Polyporales_sp"    "Peniophorella"                 
#> ASV91  NA                  "Auricularia"                   
#> ASV93  NA                  NA                              
#> ASV94  NA                  "Stereum"                       
#> ASV99  NA                  "Fomes"                         
#> ASV100 NA                  NA                              
#> ASV101 NA                  "Marchandiomyces"               
#> ASV104 NA                  "Hericium"                      
#> ASV105 NA                  "Xylodon"                       
#> ASV107 NA                  "Xylodon"                       
#> ASV108 NA                  "Exidia"                        
resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "abs_majority")@tax_table
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           ^Genus\\._abs_majority          
#> ASV7   NA                  NA                              
#> ASV8   NA                  "Stereum"                       
#> ASV12  NA                  "Xylodon"                       
#> ASV18  NA                  "Stereum"                       
#> ASV25  NA                  NA                              
#> ASV26  NA                  "Stereum"                       
#> ASV27  NA                  "Antrodiella"                   
#> ASV29  NA                  "Basidiodendron"                
#> ASV32  "Sistotrema_sp"     "Sistotrema"                    
#> ASV34  NA                  "Entocybe"                      
#> ASV35  NA                  "Fomes"                         
#> ASV41  "Mycena_sp"         "Mycena"                        
#> ASV42  NA                  "Ossicaulis"                    
#> ASV46  NA                  "Helicogloea"                   
#> ASV47  NA                  "Radulomyces"                   
#> ASV48  NA                  "Elmerina"                      
#> ASV49  NA                  "Phanerochaete"                 
#> ASV50  NA                  "Gloeohypochnicium"             
#> ASV53  NA                  "Fomes"                         
#> ASV54  NA                  "Auricularia"                   
#> ASV58  NA                  "Fomes"                         
#> ASV59  NA                  "Hyphoderma"                    
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                    
#> ASV62  "Serendipita_sp"    "Serendipita"                   
#> ASV63  NA                  NA                              
#> ASV64  NA                  "Trametes"                      
#> ASV67  NA                  "Xylodon"                       
#> ASV68  "Tricholoma_sp"     NA                              
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                  NA                              
#> ASV75  "Peniophora_reidii" "Peniophora"                    
#> ASV77  NA                  NA                              
#> ASV82  NA                  "Exidia"                        
#> ASV83  NA                  NA                              
#> ASV85  "Polyporales_sp"    NA                              
#> ASV91  NA                  "Auricularia"                   
#> ASV93  NA                  NA                              
#> ASV94  NA                  "Stereum"                       
#> ASV99  NA                  "Fomes"                         
#> ASV100 NA                  NA                              
#> ASV101 NA                  "Marchandiomyces"               
#> ASV104 NA                  "Hericium"                      
#> ASV105 NA                  "Xylodon"                       
#> ASV107 NA                  "Xylodon"                       
#> ASV108 NA                  "Exidia"                        
resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "rel_majority")@tax_table
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           ^Genus\\._rel_majority                        
#> ASV7   NA                  NA                                            
#> ASV8   NA                  "Stereum"                                     
#> ASV12  NA                  "Xylodon"                                     
#> ASV18  NA                  "Stereum"                                     
#> ASV25  NA                  "Ossicaulis/Tricholoma"                       
#> ASV26  NA                  "Stereum"                                     
#> ASV27  NA                  "Antrodiella"                                 
#> ASV29  NA                  "Basidiodendron"                              
#> ASV32  "Sistotrema_sp"     "Sistotrema"                                  
#> ASV34  NA                  "Entocybe"                                    
#> ASV35  NA                  "Fomes"                                       
#> ASV41  "Mycena_sp"         "Mycena"                                      
#> ASV42  NA                  "Ossicaulis"                                  
#> ASV46  NA                  "Helicogloea"                                 
#> ASV47  NA                  "Radulomyces"                                 
#> ASV48  NA                  "Elmerina"                                    
#> ASV49  NA                  "Phanerochaete"                               
#> ASV50  NA                  "Gloeohypochnicium"                           
#> ASV53  NA                  "Fomes"                                       
#> ASV54  NA                  "Auricularia"                                 
#> ASV58  NA                  "Fomes"                                       
#> ASV59  NA                  "Hyphoderma"                                  
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                                  
#> ASV62  "Serendipita_sp"    "Serendipita"                                 
#> ASV63  NA                  NA                                            
#> ASV64  NA                  "Trametes"                                    
#> ASV67  NA                  "Xylodon"                                     
#> ASV68  "Tricholoma_sp"     "Ossicaulis/Tricholoma"                       
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"              
#> ASV72  NA                  NA                                            
#> ASV75  "Peniophora_reidii" "Peniophora"                                  
#> ASV77  NA                  "Ossicaulis/Tricholoma"                       
#> ASV82  NA                  "Exidia"                                      
#> ASV83  NA                  NA                                            
#> ASV85  "Polyporales_sp"    "Peniophorella/Polyporales_gen_Incertae_sedis"
#> ASV91  NA                  "Auricularia"                                 
#> ASV93  NA                  NA                                            
#> ASV94  NA                  "Stereum"                                     
#> ASV99  NA                  "Fomes"                                       
#> ASV100 NA                  NA                                            
#> ASV101 NA                  "Marchandiomyces"                             
#> ASV104 NA                  "Hericium"                                    
#> ASV105 NA                  "Xylodon"                                     
#> ASV107 NA                  "Xylodon"                                     
#> ASV108 NA                  "Exidia"                                      
resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\."), method = "unanimity")@tax_table
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           ^Genus\\._unanimity             
#> ASV7   NA                  NA                              
#> ASV8   NA                  "Stereum"                       
#> ASV12  NA                  "Xylodon"                       
#> ASV18  NA                  "Stereum"                       
#> ASV25  NA                  NA                              
#> ASV26  NA                  "Stereum"                       
#> ASV27  NA                  "Antrodiella"                   
#> ASV29  NA                  "Basidiodendron"                
#> ASV32  "Sistotrema_sp"     "Sistotrema"                    
#> ASV34  NA                  "Entocybe"                      
#> ASV35  NA                  "Fomes"                         
#> ASV41  "Mycena_sp"         "Mycena"                        
#> ASV42  NA                  "Ossicaulis"                    
#> ASV46  NA                  "Helicogloea"                   
#> ASV47  NA                  "Radulomyces"                   
#> ASV48  NA                  "Elmerina"                      
#> ASV49  NA                  "Phanerochaete"                 
#> ASV50  NA                  "Gloeohypochnicium"             
#> ASV53  NA                  "Fomes"                         
#> ASV54  NA                  "Auricularia"                   
#> ASV58  NA                  "Fomes"                         
#> ASV59  NA                  "Hyphoderma"                    
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                    
#> ASV62  "Serendipita_sp"    "Serendipita"                   
#> ASV63  NA                  NA                              
#> ASV64  NA                  "Trametes"                      
#> ASV67  NA                  "Xylodon"                       
#> ASV68  "Tricholoma_sp"     NA                              
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                  NA                              
#> ASV75  "Peniophora_reidii" "Peniophora"                    
#> ASV77  NA                  NA                              
#> ASV82  NA                  "Exidia"                        
#> ASV83  NA                  NA                              
#> ASV85  "Polyporales_sp"    NA                              
#> ASV91  NA                  "Auricularia"                   
#> ASV93  NA                  NA                              
#> ASV94  NA                  "Stereum"                       
#> ASV99  NA                  "Fomes"                         
#> ASV100 NA                  NA                              
#> ASV101 NA                  "Marchandiomyces"               
#> ASV104 NA                  "Hericium"                      
#> ASV105 NA                  "Xylodon"                       
#> ASV107 NA                  "Xylodon"                       
#> ASV108 NA                  "Exidia"                        

resolve_taxo_conflict(data_fungi_mini_new, pattern_tax_ranks = c("^Genus\\.", "^Family\\.", "^Species\\."), method = "consensus")@tax_table
#> Taxonomy Table:     [45 taxa by 22 taxonomic ranks]:
#>        Domain  Phylum.x        Class.x              Order.x          
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family.x                             Genus.x             Species.x      
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                                
#> ASV7   "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV8   "Undefined Saprotroph"                                               
#> ASV12  "Undefined Saprotroph"                                               
#> ASV18  "Undefined Saprotroph"                                               
#> ASV25  "Wood Saprotroph"                                                    
#> ASV26  "Undefined Saprotroph"                                               
#> ASV27  "Wood Saprotroph"                                                    
#> ASV29  "Undefined Saprotroph"                                               
#> ASV32  "Ectomycorrhizal-Wood Saprotroph"                                    
#> ASV34  "Undefined Saprotroph"                                               
#> ASV35  "Wood Saprotroph"                                                    
#> ASV41  "Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph"
#> ASV42  "Wood Saprotroph"                                                    
#> ASV46  "Undefined Saprotroph"                                               
#> ASV47  "Undefined Saprotroph"                                               
#> ASV48  "Undefined Saprotroph"                                               
#> ASV49  "Wood Saprotroph"                                                    
#> ASV50  "Undefined Saprotroph"                                               
#> ASV53  "Wood Saprotroph"                                                    
#> ASV54  "Undefined Saprotroph"                                               
#> ASV58  "Wood Saprotroph"                                                    
#> ASV59  "Undefined Saprotroph"                                               
#> ASV61  "Undefined Saprotroph"                                               
#> ASV62  "-"                                                                  
#> ASV63  "-"                                                                  
#> ASV64  "Wood Saprotroph"                                                    
#> ASV67  "Undefined Saprotroph"                                               
#> ASV68  "Wood Saprotroph"                                                    
#> ASV71  "-"                                                                  
#> ASV72  "-"                                                                  
#> ASV75  "Plant Pathogen-Wood Saprotroph"                                     
#> ASV77  "Wood Saprotroph"                                                    
#> ASV82  "Endophyte-Undefined Saprotroph"                                     
#> ASV83  "Fungal Parasite-Undefined Saprotroph"                               
#> ASV85  "Undefined Saprotroph"                                               
#> ASV91  "Undefined Saprotroph"                                               
#> ASV93  "Wood Saprotroph-Undefined Saprotroph"                               
#> ASV94  "Undefined Saprotroph"                                               
#> ASV99  "Wood Saprotroph"                                                    
#> ASV100 "-"                                                                  
#> ASV101 "Lichen Parasite"                                                    
#> ASV104 "Undefined Saprotroph"                                               
#> ASV105 "Undefined Saprotroph"                                               
#> ASV107 "Undefined Saprotroph"                                               
#> ASV108 "Endophyte-Undefined Saprotroph"                                     
#>        Trait                 Confidence.Ranking Genus_species               
#> ASV7   "NULL"                "Probable"         "NA_NA"                     
#> ASV8   "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV12  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV18  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV25  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "White Rot"           "Probable"         "Stereum_hirsutum"          
#> ASV27  "White Rot"           "Highly Probable"  "Antrodiella_brasiliensis"  
#> ASV29  "NULL"                "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "White Rot"           "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"                "Probable"         "Entocybe_NA"               
#> ASV35  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"                "Probable"         "Mycena_renati"             
#> ASV42  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"                "Probable"         "Helicogloea_pellucida"     
#> ASV47  "White Rot"           "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"                "Probable"         "Elmerina_caryae"           
#> ASV49  "White Rot"           "Highly Probable"  "Phanerochaete_livescens"   
#> ASV50  "White Rot"           "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"                "Probable"         "Auricularia_NA"            
#> ASV58  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "White Rot"           "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "White Rot"           "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                   "-"                "NA_NA"                     
#> ASV63  "-"                   "-"                "NA_NA"                     
#> ASV64  "White Rot"           "Highly Probable"  "Trametes_versicolor"       
#> ASV67  "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV68  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                   "-"                "NA_NA"                     
#> ASV72  "-"                   "-"                "NA_NA"                     
#> ASV75  "White Rot"           "Probable"         "Peniophora_versiformis"    
#> ASV77  "Brown Rot"           "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"                "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"                "Possible"         "NA_NA"                     
#> ASV85  "White Rot"           "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"                "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"                "Probable"         "NA_NA"                     
#> ASV94  "White Rot"           "Probable"         "Stereum_ostrea"            
#> ASV99  "Brown Rot-White Rot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                   "-"                "NA_NA"                     
#> ASV101 "NULL"                "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "White Rot"           "Probable"         "Hericium_coralloides"      
#> ASV105 "White Rot"           "Probable"         "Xylodon_flaviporus"        
#> ASV107 "White Rot"           "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"                "Probable"         "Exidia_glandulosa"         
#>        Kingdom Phylum.y        Class.y           Order.y         
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"  "Cantharellales"
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV46  "Fungi" "Basidiomycota" NA                NA              
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"  "Sebacinales"   
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes" "Tremellales"   
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV83  "Fungi" "Basidiomycota" NA                NA              
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"  "Polyporales"   
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV100 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"  "Russulales"    
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"  "Agaricales"    
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"  NA              
#>        Family.y                         Genus.y                         
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  "Tricholomataceae"               "Tricholoma"                    
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  "Hydnaceae"                      "Sistotrema"                    
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  "Hyphodermataceae"               "Hyphoderma"                    
#> ASV62  "Serendipitaceae"                "Serendipita"                   
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  "Tricholomataceae"               "Tricholoma"                    
#> ASV71  "Tremellales_fam_Incertae_sedis" "Tremellales_gen_Incertae_sedis"
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  "Tricholomataceae"               "Tricholoma"                    
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  "Polyporales_fam_Incertae_sedis" "Polyporales_gen_Incertae_sedis"
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 NA                               NA                              
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species.y           ^Genus\\._consensus                           
#> ASV7   NA                  NA                                            
#> ASV8   NA                  "Stereum"                                     
#> ASV12  NA                  "Xylodon"                                     
#> ASV18  NA                  "Stereum"                                     
#> ASV25  NA                  "Ossicaulis/Tricholoma"                       
#> ASV26  NA                  "Stereum"                                     
#> ASV27  NA                  "Antrodiella"                                 
#> ASV29  NA                  "Basidiodendron"                              
#> ASV32  "Sistotrema_sp"     "Sistotrema"                                  
#> ASV34  NA                  "Entocybe"                                    
#> ASV35  NA                  "Fomes"                                       
#> ASV41  "Mycena_sp"         "Mycena"                                      
#> ASV42  NA                  "Ossicaulis"                                  
#> ASV46  NA                  "Helicogloea"                                 
#> ASV47  NA                  "Radulomyces"                                 
#> ASV48  NA                  "Elmerina"                                    
#> ASV49  NA                  "Phanerochaete"                               
#> ASV50  NA                  "Gloeohypochnicium"                           
#> ASV53  NA                  "Fomes"                                       
#> ASV54  NA                  "Auricularia"                                 
#> ASV58  NA                  "Fomes"                                       
#> ASV59  NA                  "Hyphoderma"                                  
#> ASV61  "Hyphoderma_sp"     "Hyphoderma"                                  
#> ASV62  "Serendipita_sp"    "Serendipita"                                 
#> ASV63  NA                  NA                                            
#> ASV64  NA                  "Trametes"                                    
#> ASV67  NA                  "Xylodon"                                     
#> ASV68  "Tricholoma_sp"     "Ossicaulis/Tricholoma"                       
#> ASV71  "Tremellales_sp"    "Tremellales_gen_Incertae_sedis"              
#> ASV72  NA                  NA                                            
#> ASV75  "Peniophora_reidii" "Peniophora"                                  
#> ASV77  NA                  "Ossicaulis/Tricholoma"                       
#> ASV82  NA                  "Exidia"                                      
#> ASV83  NA                  NA                                            
#> ASV85  "Polyporales_sp"    "Peniophorella/Polyporales_gen_Incertae_sedis"
#> ASV91  NA                  "Auricularia"                                 
#> ASV93  NA                  NA                                            
#> ASV94  NA                  "Stereum"                                     
#> ASV99  NA                  "Fomes"                                       
#> ASV100 NA                  NA                                            
#> ASV101 NA                  "Marchandiomyces"                             
#> ASV104 NA                  "Hericium"                                    
#> ASV105 NA                  "Xylodon"                                     
#> ASV107 NA                  "Xylodon"                                     
#> ASV108 NA                  "Exidia"                                      
#>        ^Family\\._consensus                                               
#> ASV7   "Stereaceae"                                                       
#> ASV8   "Stereaceae"                                                       
#> ASV12  "Schizoporaceae"                                                   
#> ASV18  "Stereaceae"                                                       
#> ASV25  "Lyophyllaceae/Tricholomataceae"                                   
#> ASV26  "Stereaceae"                                                       
#> ASV27  "Steccherinaceae"                                                  
#> ASV29  "Exidiaceae"                                                       
#> ASV32  "Cantharellales_fam_Incertae_sedis/Hydnaceae"                      
#> ASV34  "Entolomataceae"                                                   
#> ASV35  "Polyporaceae"                                                     
#> ASV41  "Tricholomataceae/Mycenaceae"                                      
#> ASV42  "Lyophyllaceae"                                                    
#> ASV46  "Atractiellales_fam_Incertae_sedis"                                
#> ASV47  "Pterulaceae"                                                      
#> ASV48  "Aporpiaceae"                                                      
#> ASV49  "Phanerochaetaceae"                                                
#> ASV50  "Russulales_fam_Incertae_sedis"                                    
#> ASV53  "Polyporaceae"                                                     
#> ASV54  "Auriculariaceae"                                                  
#> ASV58  "Polyporaceae"                                                     
#> ASV59  "Hyphodermataceae"                                                 
#> ASV61  "Hyphodermataceae"                                                 
#> ASV62  "Serendipitaceae"                                                  
#> ASV63  NA                                                                 
#> ASV64  "Polyporaceae"                                                     
#> ASV67  "Schizoporaceae"                                                   
#> ASV68  "Lyophyllaceae/Tricholomataceae"                                   
#> ASV71  "Tremellales_fam_Incertae_sedis"                                   
#> ASV72  NA                                                                 
#> ASV75  "Peniophoraceae"                                                   
#> ASV77  "Lyophyllaceae/Tricholomataceae"                                   
#> ASV82  "Exidiaceae"                                                       
#> ASV83  NA                                                                 
#> ASV85  "Hymenochaetales_fam_Incertae_sedis/Polyporales_fam_Incertae_sedis"
#> ASV91  "Auriculariaceae"                                                  
#> ASV93  "Stereaceae"                                                       
#> ASV94  "Stereaceae"                                                       
#> ASV99  "Polyporaceae"                                                     
#> ASV100 NA                                                                 
#> ASV101 "Corticiaceae"                                                     
#> ASV104 "Hericiaceae"                                                      
#> ASV105 "Schizoporaceae"                                                   
#> ASV107 "Schizoporaceae"                                                   
#> ASV108 "Exidiaceae"                                                       
#>        ^Species\\._consensus          
#> ASV7   NA                             
#> ASV8   "ostrea"                       
#> ASV12  "raduloides"                   
#> ASV18  "ostrea"                       
#> ASV25  "lachnopus"                    
#> ASV26  "hirsutum"                     
#> ASV27  "brasiliensis"                 
#> ASV29  "eyrei"                        
#> ASV32  "oblongisporum/Sistotrema_sp"  
#> ASV34  NA                             
#> ASV35  "fomentarius"                  
#> ASV41  "renati/Mycena_sp"             
#> ASV42  "lachnopus"                    
#> ASV46  "pellucida"                    
#> ASV47  "molaris"                      
#> ASV48  "caryae"                       
#> ASV49  "livescens"                    
#> ASV50  "analogum"                     
#> ASV53  "fomentarius"                  
#> ASV54  NA                             
#> ASV58  "fomentarius"                  
#> ASV59  "roseocremeum"                 
#> ASV61  "setigerum/Hyphoderma_sp"      
#> ASV62  "Serendipita_sp"               
#> ASV63  NA                             
#> ASV64  "versicolor"                   
#> ASV67  "raduloides"                   
#> ASV68  "lachnopus/Tricholoma_sp"      
#> ASV71  "Tremellales_sp"               
#> ASV72  NA                             
#> ASV75  "versiformis/Peniophora_reidii"
#> ASV77  "lachnopus"                    
#> ASV82  "glandulosa"                   
#> ASV83  NA                             
#> ASV85  "pubera/Polyporales_sp"        
#> ASV91  "mesenterica"                  
#> ASV93  NA                             
#> ASV94  "ostrea"                       
#> ASV99  "fomentarius"                  
#> ASV100 NA                             
#> ASV101 "buckii"                       
#> ASV104 "coralloides"                  
#> ASV105 "flaviporus"                   
#> ASV107 "raduloides"                   
#> ASV108 "glandulosa"                   
```
