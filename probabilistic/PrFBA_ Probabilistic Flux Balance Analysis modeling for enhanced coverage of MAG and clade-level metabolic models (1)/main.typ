#import "template.typ": *
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#show figure: set block(breakable: true)

#set par(
  first-line-indent: 2em,
  justify: true,
)

// Take a look at the file `template.typ` in the file panel
// to customize this template and discover how it works.
#show: project.with(
  title: "PrFBA: Probabilistic Flux Balance Analysis modeling for enhanced coverage of MAG, ASV, and clade-level metabolic models",
  authors: (
    (name: "Andrew Freiburger", email: "afreiburger@anl.gov", affiliation: "Northwestern University"),
    (name: "Filipe Liu", email: "fliu@anl.gov", affiliation: "Argonne National Laboratory"),
    (name: "Kathleen Beilsmith", email: "kbeilsmith@anl.gov", affiliation: "Argonne National Laboratory"),
    (name: "Keith Tyo", email: "k-tyo@northwestern.edu", affiliation: "Northwestern University"),
    (name: "Chris Henry", email: "chenry@anl.gov", affiliation: "Argonne National Laboratory"),
  ),
  abstract: [
    This is the abstract.
  ],
  date: "March 29, 2024",
)


= Introduction

Microbial communities are ubiquitous @Nayfach2021AMicrobiomes and serve fundamental roles as eukaryotic symbionts @Shoaie2015QuantifyingMicrobiome, @Kumar2019ModellingMicrobiome, @Shoaie2013UnderstandingModeling and pathogens @Mcguigan2015TheLung, @Jenkinson2005OralHealth, industrial assets @Tzamali2011ACommunities, @Pettit2009MixedDiscovery and antagonists @Hillman2022Top-downOperations, and ecological agents of biogeochemical cycling @Griffiths2013InsightsCommunity, @Whitman1998Prokaryotes:Majority, @Dukas2009Deep-SeaConsortia and climate regulation @Cavicchioli2019ScientistsChange, @Orphan2001Methane-consumingAnalysis, @Glass2012TraceOxide. Microbial communities are therefore vital for research fields as diverse as medicine, basic and synthetic biology and ecology @Zomorrodi2016SyntheticApplications, materials science, and civil engineering. Yet microbial communities remain under-studied @White2015AMicrobiomes, which hinders ecological research @Prosser2007TheEcology. Understanding a microbiome requires knowledge of both which biological functions are occurring within it, to reveal the ecological role of the community, and how those functions are delegated amongst its members, to predict responses of its functions to changing community structure (e.g. a new species replacing an old species) and responses of member interactions to changing environmental conditions.

The delegation of microbiome functions across its members can be elucidated if a representation or reconstruction of each member's genome is simulated in a single metabolic model @Varma1994StoichiometricW3110; however, a fundamental challenge in studying (particularly natural) microbial communities is that many of their members are unculturable and thus cannot be sequenced in isolation. Shotgun metagenome sequencing at depth, followed by assembly and binning to derive partial and often contaminated metagenome assembled genomes (MAGs) from environmental samples, is the leading method to parse microbiome members @Chivian2023. Alternative approaches such as direct annotation and long-read sequencing are impeded by short reads and no chromosomal context that create unreliable annotations and low read counts and the multitude and small abundances of members, respectively. // citation
Microbiome MAGs are nearly never complete but they are often large enough to obtain a reliable phylogenetic placement in genomic space using average nucleotide identity (ANI) approaches like GTDBtk @Chaumeil2022, which makes them unideal for directly characterizing microbiome functions despite that they are used for this purpose @Baskaran2023. Microbiome research with MAGs has several other problems, including that often only the most abundant members are able to be reconstructed into even low-quality MAGs, and deep shotgun sequencing of microbiome samples is prohibitively both costly for a single experiment (e.g. to obtain large temporal sequences of data or sample spatial variation exhaustively) and demanding of much extracted DNA. //citation

Amplicon sequencing addresses has an opposite set of strengths and weakensses as the aforementioned shutgun-sequencing-based approaches. Amplicon sequencing approaches extract a marker gene (typically 16s for microbes) from surrounding DNA and amplify it, typically using PCR approaches. This advantageously requires only infitesimal amounts of DNA and thus captures scarce species (e.g. microbes within plant or human host cells). // citation
The amplicon sequences can then be clustered into exact sequence variants (ASVs) and annotated by phylogenetically using taxonomical classifiers. This method provides no direct information about the functions present in the environment, and phylogenetic classification is typically far lower than ANI-based approaches for even a low-quality MAG, although, function can be inferred based on the functions of reference genomes in the same taxa to which the amplicon sequence has been mapped. // citation
This approach can reveal functional delegation in a micrbiome @Yilmaz2016, with both more uncertainty and more comprehensive strain mapping than MAGs analysis.

A final approach to understand functional delegation within a complex microbiome circumvents the strain-level limitations of both MAG and amplicon-based approaches and instead focuses on higher-order clades that can be more confidently supported by MAG- and ASV-based methods. This may be the desirable approach for ecosystems with unique key-stone species who general niche functions are of interest, instead of examining specific strains in the clade @McDaniel2021, @Fodelianakis2022. A clade-level study of a micrbiome pools the functions of organisms who belong to a common clade or of reference genomes within the same clade to which an ASV was assigned when inferring its function. This approach, however, can create functionally monolithic models, since a single clade can encompass numerous ecological modalities that are difficult to capture while incentivizing the proper modal expression in a given environmental condition.

All of the aforementioned methods for studying functional delegation within a microbiome require a representative model that compensates inherent uncertainties: incompleteness and contamination for MAGs; inferred functions from close reference genomes for ASVs; and underdetermined clade models. A top-down perspective with sample-level information -- community composition, chemical composition, protein expression levels, bulk fluxes, and geophysical characteristics -- may alleviate some uncertainties that arise from bottom-up functional reconstruction. //citation
Some uncertainties can also be resolved by studying functions as working parts within a broader system. // what is meant by this?
Flux balance analysis @Orth2010, @Lee2006FluxMetabolomics and community metabolic modeling @Chan2017SteadyCom:Stability, with new formalisms that we developed to address the aforementioned uncertainties, can achieve this for metabolism.

Herein, we propose a new probabilistic approach that captures microbiome systems from MAGs, ASVs, or clades in single genome-scale metabolic models (GEMs) @Edwards2000TheCapabilities, @Collins2003TheBiology, @Covert2001MetabolicSilico, called Probabilistic Flux Balance Analysis (prFBA). Several methods have applied probability in reconstructing GEMs @Qi2014, predicting microbiome structure from chemical data @Eveillard2019, or merging regulatory networks into GEMs @Yu2019; however, there is no method that captures uncertainty in a MAG or ASV, or diversity of a phylogenetic clade, in a GEM and simulates metabolism according to these uncertainties or niche diversity.  prFBA favors the most conserved functions that are determined by the frequency, or probability, of each function among the closest reference genomes or genomes of a specified clade.  prFBA and the model for constructing probabilistic models are currently available as a ModelSEED Python API for local developer access, but will be expanded for boarder usage through KBase Applications. The probabilistic modeling framework presented here -- both the construction of probabilistic models and the prFBA simulation method -- will accelerate microbial ecology and offers a unique method for accommodating uncertainty in MAGs and ASV, and uniquely facilitate clade-level microbiome studies.

// At the same time, the next-gen sequencing revolution has yielded hundreds of thousands of reference genomes that are complete and high quality. We can leverage these reference genomes to fill gaps in our MAGs so long as we also carefully acknowledge and address the uncertainty that arises in doing so. 

// NOTE FROM CHRIS: I suggest cutting the next three sections, merging them when appropriate with my text above. Generally introductions do not have subsections like this. I suggest adding you probabilistic modeling content to the end of the text above and smoothing things out.

// cite former preprints/papers where we elaborate the essentiality of communities

// MAGs are challenging because they are fragmentary and low-quality, so anchoring them to reference genomes is essential to acquire the most complete assessment.

// mine Dylan's paper for citations and content pertaining to MAGs



= Methods 

== Mapping genome or MAG in phylogenetic space <mapping_method>

Our method begins with the parameterization of a genome or MAG sequence, which is localized within phylogenetic space of a reference database or MAG collection. We used the Genome Taxonomical Database (GTDB) @Parks2022 because it covers vast phylogenetic space and has a convenient tool-kit API @Chaumeil2022 that is integrated into KBase for rapid and reproducible computations @Arkin2018. The threshold for acquiring the set of is //

== Constructing a probabilistic model <construction>

The creation of probabilistic models from a collection of organisms is performed by determining the frequency of a given model attribute: annotation, gene, or reaction, in order of preference.  The mapping of annotations or genes is ideal because it can occur before the reconstruction of models into GEMs, so only the probabilistic model is reconstructed from the union of all annotations or genes with the respective probabilities. The mapping of reactions, by contrast, requires each organism to be reconstructed and then the probabilistic model is created from the union of organism reactions with their respective frequencues. The probability for a given function object

$ p_"func" = "num"_"organisms,func" / "num"_"organisms" $ <frequency>

is determined as its frequency of the among the set of reference organisms, and is stored as an attribute of each respective reaction. The probability defined in *@frequency* can be further refined for ASV and MAG systems by applying a secondary weighting $alpha$ that favor organisms in the pool of MAG or ASV organisms who are closer to the reference organism. The $alpha_"i"$ of organism $i$

$ alpha_"i" = ("(" "num"_"point,func" #sym.sect "num"_"i,func" ")" * "completeness") / "num"_"i,func" $

is the fraction of functions that are shared between the MAG or ASV point and the phylogenetically close reference organism $i$, according to ANI distance, multiplied by the completeness of the MAG or ASV. The $0<alpha<1$ is then used to compute the $0<p_"func"<1$ function probability

$ p_"func" = (sum_i^I (b_"func,i" * alpha_"i") )/ (sum_"i"^"I" (alpha_"i")) $

for all genomes $i$ that are close to the examined $"org"$, where $b_"func,i"$ is a binary variable that indicates whether the function is (1) or is not (0) in the examined genome $i$. 

// The updated determination of probability

// $ p_"func" = ("num"_"organisms,func"*alpha_"organism") / "num"_"organisms" $

// thereby approximates organism behavior based on the closest reference genomes.

// probabilistic annotations are the ideal method for constructing a mode, e.g. through the MS pipeline, but can also be applied for existing models, such as AGORA2, via the reactions.

== Simulating probabilistic FBA (prFBA)

The prFBA method is defined by several constraints, simulations, and a unique objective function. One of the constraints ("ElementUptake") limits carbon consumption to a specified $"ele"_"limit"$ number of atoms

$ sum^"EX"_"ex" ("ele"_"ex"*("ex"_"forwardVar" #sym.plus.circle "ex"_"reverseVar")) = "ele"_"limit" $

which imposes efficient utilization of resources from a given media, where complete media is the default environment. This constraint requires the atom count for each of the exchanged metabolites $"ele"_"ex"$, which is then multipled by either the forward or reverse flux variables $"ex"_"forwardVar" #sym.plus.circle "ex"_"ReverseVar"$ when calculating the impact of the variable on the total element balance. A second constraint ("CommKinetics"), limits the total reaction flux of each model to a $"kinCoef"$ multiple of its biomass growth

$ sum^"R"_"r" ("r"_"forwardVar"+"r"_"reverseVar") = "kinCoef" * "bio"_"forwardVar"  $

which effectively prevents organisms from being metabolically exploited for efficient reactions or pathways without themselves growing. The final constraint ("min_biomass") simulates the model with the default objective function of maximizing biomass growth and mandates that the biomass flux is at least 95% of its maximum value

$ "cons": "bio1"_"flux" > 0.95*"bio1"_"flux, max" $

to maintain optimum cellular growth despite later replacing the objective function.



The prFBA objective function

$ min{r,"ex"}: sum_r^"R"_"internal" ((1-p^"prob_exp")*r_"flux")+"min_prob") + \ sum_"rd"^"RD" (-"rd"_"expression"*"rd"_"flux")) + sum_"ex"^"EX" (100*"ex"_"flux") $

finds the prototypical metabolic activity of a MAG, clade, or ASV through three sums. The first sum penalizes intracellular reactions, for which data does not exist, inversely with their probability by scaling their flux $r_"flux"$. The $"min_prob"$ parameter establishes the smallest probability that is assigned to a reaction, which prevents some reactions from becoming penalized into oblivion. The $"prob_exp"$ parameter tailors the significance $p_"object"$ in determining the reaction weightings. The second sum preferentially rewards internal reactions for which expression data exists ($"rd"$) proportionally with the relative expression of each respective reaction ($"rd"_"expression"$). The third sum penalizes exchanges, which maximizes transport and thus interactions among community members. This objective further mitigates degenerate solutions because alternative reaction pathways will nearly always have different weightings and thus the minimization will just use the most probable pathway.



== Data integration

Meta-transcriptomics data can be invaluable for this method by yielding a probability of 1 for genes that are expressed in the data, because the uncertainty in gene expression that is captured by $p_"object"$ disappears.



= Results

== New Probabilistic Flux Balance Analysis Approach to GEM Reconstruction and Analysis

Our workflow illustrated in *@workflow* and analyzes GEMs for individual strains with an imcomplete sequence (MAG or ASV) or taxa with exemplar genomes, as a compliment to ModelSEED2 @Faria2023 or CarveME @Machado2018 and COBRApy @Ebrahim2013 that build and analyze GEMs from complete high-quality genome sequences, respectively. Our approach reconstructs strain or clade probabilistic models by mapping a MAG, ASV, or taxonomical class into phylogenetic space (see *@mapping_method*). // this needs to be created, by Chris?
The phylogenetically closest reference genomes or other MAGs are then identified, preferably being with genomes of >0.95 ANI score but ultimately including the closest representative genomes because MAGs and ASVs may not have very close reference genomes. Fortunately, our proposed method is somewhat flexible to the distance of the reference genomes, where the uncertainty is likely greater for each of the models.

All nearby genomes are then annotated with a model-friendly functional ontology such as RAST, GO, and KEGGKO.  We employ RAST because it is compatible with our ModelSEED2 reconstruction pipeline and because we did not focus on annotation uncertainty in this work, but we strongly recommend aggregating the genome annotations from through several methods since functional annotation significantly contributes to uncertainty in these models. The annotated genomes of the selected genomes are then merged into a consolidated probabilitic annotation (CPA, see *@construction*) that includes all unique functions among the selected genomes and is each assigned a probability based on two factors: (1) how prevalent the function is amongst the selected genomes included in our set; and (2) how distant each genome containing the function is from the MAGs and ASV sequence. 

Our approach also applies to existing metabolic models of reference genomes by merging their biochemical reactions as a proxy of functions instead of the original annotations, which is demonstrated by representing a gut microbiome with the AGORA2 model collection. 

A reconstruction method is next applied based on the CPA, or the CPA becomes the GEM itself when merging existing metabolic models. We used the ModelSEED2 method to build our CPA model for the MAG and clade systems, and ascribed our function-level probabilities to all GEM reactions. CPA-constructed GEMS can be simulated individually or aggregated into a compartmentalized community metabolic model to probabilistically study a microbiome system where each probabilistic model is constrained based on species abundance of the microbiome sample. Our example here all study microbiomes through CPA-constructed community metabolic models as the more complex use-case of probabilistic modeling.

Finally, flux balance analysis with numerous potential constraints and modeling modalities is employed to predict metabolic activity and fluxes that maximally use the highest probability model reactions, which captures as much top-down knowledge of the system behavior such as gene/protein expression, species abundance, global flux measures, metabolomics profiles, QSIP profiles, and TN-seq data. Low probability reactions will generally be avoided, but could still be used if they are essential for explaining systems-level experimental observations, which allows CPA-constructed models, compartmentalized community models of these individual models, are able to resolve uncertainties based on systems level observations.

= Conclusion

Probabilistic GEM modeling, introduced herein, should be powerful tool for microbial ecology. A MAG from an unstudied soil sample, for example, could map to the Psueodmonad genus but an insufficient amount of genetic sequence is available to map to a species. Our method would be able to assemble a genome-scale metabolic model that pulls reactions from the closest reference genomes and would exercise the most conserved behavior of the phylogenetic space while implementing possible niche behavior as environmental constraints necessitate.

Our method further uniquely permits an analysis of prototypical clade behavior and specifically how multiple clades interact with each other. This introduces a new paradigm where the general roles of microbes can be contextualized in otherwise excessively complex ecological systems. 


#columns(2, gutter:10pt)[
  #figure(
    image("methods_figure.jpg"),
    caption: [
      Panel a denotes the genome space defined by reference genomes and either a MAG, ASV, or clade of interest.  The MAG and ASV sequences are mapped into this space based on average nucleotide identity (ANI). The phylogenetically close reference genomes and (other) MAGs can then be used for gapfilling and model reconstruction by merging the annotations. The development of a clade model simply assigns a cut-off at a taxonomical level and involves merging the reactions per se from the set of all models.
    ], // the caption is too big to fit in the bottom margin of the figure, and thus requires wrapping over
       // two columns and outside of the caption attribute in the below paragraphs.
  ) <workflow>
  Panel b depicts three reference organisms that are either phylogenetically close to a MAG/ASV or are in a common clade, and are merged into a single model that captures metabolic uncertainty and diversity from the closest reference organisms by assigning a probability to each reaction according to its frequency among the reference organisms. Reactions that are already annotated to the MAG or are empirically detected are overwritten with a probability of 1. The red sections of the model Venn diagram are the reactions that are present in only one of the three models; the green sections are the reactions that are present in two of the three models; and the teal center are the reactions present in all of the models.
      
  Panel c illustrates a sample PrFBA simulation of a probabilistic model.   Example alternative pathways from the 3 reference models that convert compound 1 (cpd1) into compound 5 (cpd5) in the left columns are examined probabilistically in the right diagram, where the reaction arrows are painted with the respective probabilities of each reaction according to the color schema from panel b. PrFBA would preferentially follow the pathway $"cpd1" => "cpd2" => "cpd3" => "cpd4" => "cpd5"$ because this is the route that employs the most probable reactions, unless additional constraints incentivize the use of less probable reactions. This probabilistic framework allows a broad assessment of MAG, ASV, or clade behavior in a single-strain GEM. The probabilistic model can also model sample-level systems, where all possible hypotheses of functional content in a MAG are evaluated to discern how biological activity is delegated across a microbial community.
]


// The probabilistic model and PrFBA methods are successfully integrated into KBase Applications, exhibited in *Figure 2*.  The probabilistic model construction App accepts genomes or a genomeset and assigns probabilities to functions at the genome-level, instead of the reactions-level of completed FBA reactions, to expedite the process.  The probabilistic set of functions are then reconstructed into a metabolic model with stored probabilities in the reaction attributes.  

// == KBase Applications


// == Scientific insights

// what is the interpretation of the results from Kayla's simulations



// provide the bibliograpy of the 
#bibliography("probabilistic.bib")
