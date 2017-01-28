# Thermodynamic Promoter Activity Predictor (TPAP)

## Background
Most DNA-protein binding prediction models for promoter association tend to work by comparing amino acid sequence and tertiary structure homology in order to find potential transcription factors and binding sites.
However, these prediction models cannot tell us anything about the overall change in affinity between conserved promoter regions and transcription factors,
and its consequential effect on promoter activity when point mutations are introduced into the promoter. Creating complex prediction models for this type of 
assesment requires huge experimental effort (e.g. numerous promoter-GFP runs with high sequence and protein complex specificity), and so computational models are needed to diminish this work load
and to better understand promoter activity for future applications. This project's purpose is to provide a framework for experimental research groups that enables them to train a model using experimental data capable
of determining effects in promoter activity as inserted mutations alter specific and/or non-specific binding sites in a promoter.
Although this project is at its infancy, the more experimental data sets provided to the model, the more accurate future predictions will be - suggesting that the growth and availability of projects like this are heavily influenced by communal research. 

## Theory
There are many ways to try and understand gene regulation networks, but it can be said that all regulation stages are determined initially by the
probability of a determined number of molecules to be bound to their respective target regions. For RNA Polymerase-Promoter association, the statistical approximation of binding probability can be represented as follows: 

![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/8ff0e6ab5216ce921397818c521c2631a12bd50f) 

where  ![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/b4dc73bf40314945ff376bd363916a738548d40a) = effective number of RNA Polymerase molecules available for binding (tricky!),

 ![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/c6112a38cefce408105ca77bf32a8df4898de065) = non-specific sites available (i.e. bp not participating in gene expression, tricky in eukaryotes),

![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/6a921feebd7e1d5427126b3b7d1d68e5313fc7e2) (E<sub>S</sub> = Mean binding energy of specific sites, E<sub>NS</sub> = Mean binding energy of non-specific sites),

and ![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/36c396599916ca276daed739acd8df48718641cb) = Average impact of associated transcriptional factors.


From this it becomes quite obvious that the most important factor in this mathematical approach is determining a 'good enough' value for ![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/36c396599916ca276daed739acd8df48718641cb).
But how can we measure a value for ![equation](https://wikimedia.org/api/rest_v1/media/math/render/svg/36c396599916ca276daed739acd8df48718641cb), and when (or how) do we deem values as 'good enough' (i.e. empirically acceptable)?
One approach used in the past to assay increase or decrease in gene expression based on TF binding consists in the use of Mixed Interger Linear Programming (MILP) models that produce Steiner trees of biochemical signalling networks in order to find transcriptional patterns of interest.
Combining this approach with a learning model that considers previously observed effects in target genes (i.e. 1 for increase in transcription, 0 for decrease in transcription) offers a rough estimate of the overall effects in promoter activity due to succesful individual promoter-TF binding events.  

## Databases

## Author
- [Charles Sanfiorenzo Cruz]() - charles.sanfiorenzo@upr.edu

## References
-  Bintu, L.; Buchler, N; Garcia, H; Gerland, U et al. (2005). "Transcriptional regulation by the numbers: models". Current Opinion in Genetics & Development. 15: 116–124
