# Compute congruence metrics between two taxonomic assignments

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Computes metrics quantifying the congruence between two taxonomic
assignments for the same set of taxa (ASVs/OTUs). This is useful for
comparing different taxonomic databases, assignment methods, or
reference versions.

## Usage

``` r
tc_congruence_metrics(physeq_1, physeq_2 = NULL, ranks_1, ranks_2 = NULL)
```

## Arguments

- physeq_1:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with the first taxonomic assignment.

- physeq_2:

  (phyloseq, default NULL) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with the second taxonomic assignment. If NULL, uses physeq_1
  (useful for comparing different rank columns from the same object).

- ranks_1:

  (character vector, required) Taxonomic rank names to use for the first
  assignment. Must match column names in physeq_1's tax_table.

- ranks_2:

  (character vector, default NULL) Taxonomic rank names to use for the
  second assignment. Must match column names in physeq_2's tax_table. If
  NULL, uses the same ranks as ranks_1 vector.

## Value

A list with the following components:

- summary:

  A data frame with counts and percentages for each category, including
  `leaf_match_congruence` which counts taxa where the deepest
  classification matches regardless of the path

- only_in_1:

  Character vector of taxa names present only in physeq_1

- only_in_2:

  Character vector of taxa names present only in physeq_2

- classified_only_1:

  Character vector of common taxa classified only in physeq_1 (all NA in
  physeq_2 for the specified ranks)

- classified_only_2:

  Character vector of common taxa classified only in physeq_2 (all NA in
  physeq_1 for the specified ranks)

- unclassified_both:

  Character vector of common taxa with all NA in both physeq objects

- total_congruent:

  Character vector of taxa with identical taxonomic paths (same values
  at all ranks, same NA positions)

- partial_1_deeper:

  Character vector of taxa where physeq_1 classifies to a deeper rank
  but agrees on shared ranks

- partial_2_deeper:

  Character vector of taxa where physeq_2 classifies to a deeper rank
  but agrees on shared ranks

- incongruent_leaves:

  Character vector of taxa with disagreement at the deepest (leaf) level
  where both have assignments

- incongruent_nodes:

  Character vector of taxa with disagreement at higher (internal node)
  levels, even if leaves match

- details:

  A data frame with per-taxon details including: taxon name, depth_1,
  depth_2, leaf_1, leaf_2, leaf_match (TRUE if leaves are equal), and
  category

## Author

Adrien Taudière

## Examples

``` r
# Compare two taxonomic assignments from the same phyloseq object
metrics <- tc_congruence_metrics(
  subset_taxa(Glom_otu, Phyla == "Fungi" | Phylum__eukaryome_Glomero == "Fungi"),
  ranks_1 = c("Phyla", "Class", "Order", "Family"),
  ranks_2 = c(
    "Phylum__eukaryome_Glomero", "Class__eukaryome_Glomero",
    "Order__eukaryome_Glomero", "Family__eukaryome_Glomero"
  )
)

# View summary
metrics$summary
#>                 category count percentage
#> 1              only_in_1     0       0.00
#> 2              only_in_2     0       0.00
#> 3      classified_only_1     0       0.00
#> 4      classified_only_2     0       0.00
#> 5      unclassified_both     0       0.00
#> 6        total_congruent     0       0.00
#> 7  leaf_match_congruence   111       9.68
#> 8       partial_1_deeper     0       0.00
#> 9       partial_2_deeper     0       0.00
#> 10    incongruent_leaves   673      58.67
#> 11     incongruent_nodes   474      41.33

# Get taxa with leaf-level incongruence
metrics$incongruent_leaves
#>   [1] "Taxa_1"    "Taxa_8"    "Taxa_9"    "Taxa_14"   "Taxa_15"   "Taxa_18"  
#>   [7] "Taxa_19"   "Taxa_27"   "Taxa_29"   "Taxa_30"   "Taxa_31"   "Taxa_35"  
#>  [13] "Taxa_37"   "Taxa_52"   "Taxa_57"   "Taxa_58"   "Taxa_60"   "Taxa_61"  
#>  [19] "Taxa_65"   "Taxa_66"   "Taxa_67"   "Taxa_68"   "Taxa_69"   "Taxa_70"  
#>  [25] "Taxa_72"   "Taxa_73"   "Taxa_77"   "Taxa_79"   "Taxa_80"   "Taxa_86"  
#>  [31] "Taxa_87"   "Taxa_90"   "Taxa_94"   "Taxa_95"   "Taxa_101"  "Taxa_102" 
#>  [37] "Taxa_106"  "Taxa_108"  "Taxa_109"  "Taxa_110"  "Taxa_111"  "Taxa_113" 
#>  [43] "Taxa_115"  "Taxa_123"  "Taxa_126"  "Taxa_128"  "Taxa_129"  "Taxa_132" 
#>  [49] "Taxa_136"  "Taxa_137"  "Taxa_139"  "Taxa_140"  "Taxa_143"  "Taxa_145" 
#>  [55] "Taxa_147"  "Taxa_148"  "Taxa_151"  "Taxa_157"  "Taxa_166"  "Taxa_175" 
#>  [61] "Taxa_179"  "Taxa_181"  "Taxa_184"  "Taxa_185"  "Taxa_187"  "Taxa_194" 
#>  [67] "Taxa_195"  "Taxa_200"  "Taxa_204"  "Taxa_207"  "Taxa_214"  "Taxa_222" 
#>  [73] "Taxa_223"  "Taxa_228"  "Taxa_229"  "Taxa_234"  "Taxa_235"  "Taxa_236" 
#>  [79] "Taxa_246"  "Taxa_248"  "Taxa_252"  "Taxa_257"  "Taxa_259"  "Taxa_263" 
#>  [85] "Taxa_266"  "Taxa_269"  "Taxa_271"  "Taxa_274"  "Taxa_285"  "Taxa_288" 
#>  [91] "Taxa_289"  "Taxa_292"  "Taxa_293"  "Taxa_294"  "Taxa_301"  "Taxa_304" 
#>  [97] "Taxa_305"  "Taxa_309"  "Taxa_310"  "Taxa_312"  "Taxa_313"  "Taxa_316" 
#> [103] "Taxa_317"  "Taxa_318"  "Taxa_322"  "Taxa_328"  "Taxa_329"  "Taxa_339" 
#> [109] "Taxa_341"  "Taxa_342"  "Taxa_345"  "Taxa_346"  "Taxa_348"  "Taxa_352" 
#> [115] "Taxa_354"  "Taxa_357"  "Taxa_358"  "Taxa_362"  "Taxa_364"  "Taxa_367" 
#> [121] "Taxa_378"  "Taxa_379"  "Taxa_390"  "Taxa_402"  "Taxa_403"  "Taxa_407" 
#> [127] "Taxa_408"  "Taxa_412"  "Taxa_418"  "Taxa_419"  "Taxa_420"  "Taxa_422" 
#> [133] "Taxa_427"  "Taxa_433"  "Taxa_440"  "Taxa_443"  "Taxa_450"  "Taxa_453" 
#> [139] "Taxa_454"  "Taxa_459"  "Taxa_460"  "Taxa_465"  "Taxa_480"  "Taxa_483" 
#> [145] "Taxa_484"  "Taxa_485"  "Taxa_489"  "Taxa_498"  "Taxa_503"  "Taxa_507" 
#> [151] "Taxa_509"  "Taxa_510"  "Taxa_511"  "Taxa_521"  "Taxa_526"  "Taxa_532" 
#> [157] "Taxa_535"  "Taxa_541"  "Taxa_547"  "Taxa_549"  "Taxa_551"  "Taxa_556" 
#> [163] "Taxa_557"  "Taxa_558"  "Taxa_561"  "Taxa_569"  "Taxa_570"  "Taxa_574" 
#> [169] "Taxa_583"  "Taxa_585"  "Taxa_588"  "Taxa_596"  "Taxa_597"  "Taxa_598" 
#> [175] "Taxa_599"  "Taxa_600"  "Taxa_603"  "Taxa_604"  "Taxa_605"  "Taxa_608" 
#> [181] "Taxa_610"  "Taxa_611"  "Taxa_619"  "Taxa_625"  "Taxa_626"  "Taxa_630" 
#> [187] "Taxa_632"  "Taxa_635"  "Taxa_640"  "Taxa_648"  "Taxa_650"  "Taxa_658" 
#> [193] "Taxa_659"  "Taxa_660"  "Taxa_662"  "Taxa_665"  "Taxa_668"  "Taxa_669" 
#> [199] "Taxa_670"  "Taxa_673"  "Taxa_682"  "Taxa_695"  "Taxa_696"  "Taxa_702" 
#> [205] "Taxa_707"  "Taxa_708"  "Taxa_709"  "Taxa_711"  "Taxa_715"  "Taxa_729" 
#> [211] "Taxa_730"  "Taxa_733"  "Taxa_745"  "Taxa_747"  "Taxa_751"  "Taxa_755" 
#> [217] "Taxa_759"  "Taxa_766"  "Taxa_769"  "Taxa_774"  "Taxa_776"  "Taxa_789" 
#> [223] "Taxa_795"  "Taxa_796"  "Taxa_797"  "Taxa_806"  "Taxa_807"  "Taxa_809" 
#> [229] "Taxa_811"  "Taxa_820"  "Taxa_826"  "Taxa_828"  "Taxa_850"  "Taxa_856" 
#> [235] "Taxa_860"  "Taxa_865"  "Taxa_866"  "Taxa_876"  "Taxa_881"  "Taxa_888" 
#> [241] "Taxa_898"  "Taxa_899"  "Taxa_901"  "Taxa_902"  "Taxa_903"  "Taxa_932" 
#> [247] "Taxa_933"  "Taxa_934"  "Taxa_950"  "Taxa_955"  "Taxa_966"  "Taxa_968" 
#> [253] "Taxa_969"  "Taxa_971"  "Taxa_973"  "Taxa_983"  "Taxa_992"  "Taxa_995" 
#> [259] "Taxa_996"  "Taxa_1001" "Taxa_1008" "Taxa_1020" "Taxa_1021" "Taxa_1023"
#> [265] "Taxa_1026" "Taxa_1031" "Taxa_1033" "Taxa_1036" "Taxa_1041" "Taxa_1055"
#> [271] "Taxa_1071" "Taxa_1074" "Taxa_1075" "Taxa_1086" "Taxa_1093" "Taxa_1101"
#> [277] "Taxa_1113" "Taxa_1121" "Taxa_1125" "Taxa_1135" "Taxa_1142" "Taxa_1150"
#> [283] "Taxa_1152" "Taxa_1160" "Taxa_1161" "Taxa_1167" "Taxa_1168" "Taxa_1180"
#> [289] "Taxa_1187" "Taxa_1193" "Taxa_1195" "Taxa_1210" "Taxa_1217" "Taxa_1233"
#> [295] "Taxa_1236" "Taxa_1248" "Taxa_1252" "Taxa_1258" "Taxa_1263" "Taxa_1265"
#> [301] "Taxa_1274" "Taxa_1276" "Taxa_1277" "Taxa_1279" "Taxa_1286" "Taxa_1294"
#> [307] "Taxa_1296" "Taxa_1297" "Taxa_1304" "Taxa_1311" "Taxa_1327" "Taxa_1329"
#> [313] "Taxa_1335" "Taxa_1353" "Taxa_1374" "Taxa_1381" "Taxa_1382" "Taxa_1385"
#> [319] "Taxa_1387" "Taxa_1388" "Taxa_1390" "Taxa_1391" "Taxa_1424" "Taxa_1430"
#> [325] "Taxa_1432" "Taxa_1433" "Taxa_1436" "Taxa_1437" "Taxa_1449" "Taxa_1453"
#> [331] "Taxa_1465" "Taxa_1471" "Taxa_1474" "Taxa_1481" "Taxa_1484" "Taxa_1499"
#> [337] "Taxa_1501" "Taxa_1504" "Taxa_1521" "Taxa_1527" "Taxa_1529" "Taxa_1536"
#> [343] "Taxa_1546" "Taxa_1553" "Taxa_1558" "Taxa_1564" "Taxa_1573" "Taxa_1585"
#> [349] "Taxa_1591" "Taxa_1600" "Taxa_1608" "Taxa_1629" "Taxa_1635" "Taxa_1660"
#> [355] "Taxa_1681" "Taxa_1693" "Taxa_1699" "Taxa_1711" "Taxa_1719" "Taxa_1723"
#> [361] "Taxa_1741" "Taxa_1745" "Taxa_1756" "Taxa_1768" "Taxa_1770" "Taxa_1783"
#> [367] "Taxa_1793" "Taxa_1795" "Taxa_1804" "Taxa_1818" "Taxa_1821" "Taxa_1823"
#> [373] "Taxa_1830" "Taxa_1832" "Taxa_1834" "Taxa_1852" "Taxa_1853" "Taxa_1857"
#> [379] "Taxa_1873" "Taxa_1874" "Taxa_1879" "Taxa_1882" "Taxa_1892" "Taxa_1898"
#> [385] "Taxa_1904" "Taxa_1909" "Taxa_1912" "Taxa_1914" "Taxa_1915" "Taxa_1920"
#> [391] "Taxa_1930" "Taxa_1960" "Taxa_1962" "Taxa_1977" "Taxa_1988" "Taxa_1999"
#> [397] "Taxa_2008" "Taxa_2021" "Taxa_2044" "Taxa_2047" "Taxa_2054" "Taxa_2072"
#> [403] "Taxa_2081" "Taxa_2102" "Taxa_2108" "Taxa_2114" "Taxa_2118" "Taxa_2131"
#> [409] "Taxa_2134" "Taxa_2143" "Taxa_2155" "Taxa_2164" "Taxa_2171" "Taxa_2172"
#> [415] "Taxa_2201" "Taxa_2205" "Taxa_2227" "Taxa_2239" "Taxa_2242" "Taxa_2247"
#> [421] "Taxa_2253" "Taxa_2256" "Taxa_2270" "Taxa_2279" "Taxa_2297" "Taxa_2299"
#> [427] "Taxa_2302" "Taxa_2304" "Taxa_2315" "Taxa_2316" "Taxa_2342" "Taxa_2351"
#> [433] "Taxa_2360" "Taxa_2367" "Taxa_2371" "Taxa_2374" "Taxa_2386" "Taxa_2411"
#> [439] "Taxa_2415" "Taxa_2426" "Taxa_2435" "Taxa_2440" "Taxa_2443" "Taxa_2444"
#> [445] "Taxa_2472" "Taxa_2493" "Taxa_2530" "Taxa_2536" "Taxa_2541" "Taxa_2549"
#> [451] "Taxa_2553" "Taxa_2559" "Taxa_2566" "Taxa_2574" "Taxa_2589" "Taxa_2603"
#> [457] "Taxa_2611" "Taxa_2614" "Taxa_2626" "Taxa_2631" "Taxa_2637" "Taxa_2639"
#> [463] "Taxa_2640" "Taxa_2644" "Taxa_2654" "Taxa_2663" "Taxa_2680" "Taxa_2682"
#> [469] "Taxa_2698" "Taxa_2706" "Taxa_2708" "Taxa_2734" "Taxa_2742" "Taxa_2749"
#> [475] "Taxa_2754" "Taxa_2766" "Taxa_2768" "Taxa_2777" "Taxa_2778" "Taxa_2780"
#> [481] "Taxa_2793" "Taxa_2807" "Taxa_2811" "Taxa_2812" "Taxa_2813" "Taxa_2818"
#> [487] "Taxa_2819" "Taxa_2830" "Taxa_2833" "Taxa_2834" "Taxa_2839" "Taxa_3201"
#> [493] "Taxa_2857" "Taxa_2859" "Taxa_2861" "Taxa_2870" "Taxa_2875" "Taxa_2885"
#> [499] "Taxa_2886" "Taxa_2889" "Taxa_2902" "Taxa_2904" "Taxa_2906" "Taxa_2908"
#> [505] "Taxa_2911" "Taxa_2912" "Taxa_2929" "Taxa_2930" "Taxa_2934" "Taxa_2938"
#> [511] "Taxa_2945" "Taxa_2951" "Taxa_2954" "Taxa_2962" "Taxa_2970" "Taxa_2973"
#> [517] "Taxa_2998" "Taxa_2999" "Taxa_3000" "Taxa_3021" "Taxa_3024" "Taxa_3031"
#> [523] "Taxa_3041" "Taxa_3045" "Taxa_3047" "Taxa_3048" "Taxa_3055" "Taxa_3059"
#> [529] "Taxa_3073" "Taxa_3087" "Taxa_3099" "Taxa_3103" "Taxa_3108" "Taxa_3110"
#> [535] "Taxa_3118" "Taxa_3124" "Taxa_3134" "Taxa_3141" "Taxa_3149" "Taxa_3157"
#> [541] "Taxa_3158" "Taxa_3188" "Taxa_3190" "Taxa_3198" "Taxa_3207" "Taxa_3209"
#> [547] "Taxa_3213" "Taxa_3216" "Taxa_3217" "Taxa_3220" "Taxa_3229" "Taxa_3230"
#> [553] "Taxa_3244" "Taxa_3264" "Taxa_3265" "Taxa_3267" "Taxa_3274" "Taxa_3275"
#> [559] "Taxa_3276" "Taxa_3277" "Taxa_3278" "Taxa_3290" "Taxa_3291" "Taxa_3300"
#> [565] "Taxa_3308" "Taxa_3310" "Taxa_3312" "Taxa_3314" "Taxa_3319" "Taxa_3328"
#> [571] "Taxa_3331" "Taxa_3338" "Taxa_3339" "Taxa_3341" "Taxa_3342" "Taxa_3346"
#> [577] "Taxa_3347" "Taxa_3348" "Taxa_3355" "Taxa_3356" "Taxa_3359" "Taxa_3360"
#> [583] "Taxa_3361" "Taxa_3363" "Taxa_3365" "Taxa_3372" "Taxa_3375" "Taxa_3379"
#> [589] "Taxa_3387" "Taxa_3388" "Taxa_3389" "Taxa_3390" "Taxa_3394" "Taxa_3398"
#> [595] "Taxa_3401" "Taxa_3402" "Taxa_3403" "Taxa_3405" "Taxa_3415" "Taxa_3417"
#> [601] "Taxa_3422" "Taxa_3428" "Taxa_3436" "Taxa_3438" "Taxa_3439" "Taxa_3445"
#> [607] "Taxa_3446" "Taxa_3447" "Taxa_3450" "Taxa_3451" "Taxa_3455" "Taxa_3456"
#> [613] "Taxa_3458" "Taxa_3465" "Taxa_3466" "Taxa_3469" "Taxa_3470" "Taxa_3471"
#> [619] "Taxa_3474" "Taxa_3478" "Taxa_3482" "Taxa_3486" "Taxa_3488" "Taxa_3490"
#> [625] "Taxa_3493" "Taxa_3496" "Taxa_3506" "Taxa_3509" "Taxa_3512" "Taxa_3529"
#> [631] "Taxa_3531" "Taxa_3549" "Taxa_3550" "Taxa_3552" "Taxa_3558" "Taxa_3566"
#> [637] "Taxa_3567" "Taxa_3571" "Taxa_3572" "Taxa_3574" "Taxa_3576" "Taxa_3577"
#> [643] "Taxa_3578" "Taxa_3602" "Taxa_3613" "Taxa_3614" "Taxa_3619" "Taxa_3620"
#> [649] "Taxa_3622" "Taxa_3624" "Taxa_3626" "Taxa_3629" "Taxa_3635" "Taxa_3637"
#> [655] "Taxa_3640" "Taxa_3645" "Taxa_3646" "Taxa_3647" "Taxa_3649" "Taxa_3650"
#> [661] "Taxa_3654" "Taxa_3661" "Taxa_3664" "Taxa_3666" "Taxa_3671" "Taxa_3672"
#> [667] "Taxa_3673" "Taxa_3676" "Taxa_3677" "Taxa_3678" "Taxa_3692" "Taxa_3693"
#> [673] "Taxa_3697"
```
