#Now building clade genome objects and saving to file (will be loaded to KBase)
from cobrakbase.core.kbasegenome.genome import KBaseGenome

mags = util.load("mag_list_v2")
functions = {}
features = {}
feature_aliases = {}
feature_probabilities = {}
md5_list = []
count = 0
for mag in mags:
    #Loading template genomex
    mag = list(map(str, mag))
    # print(mag)
    genome = util.get_object('/'.join([mag[6], mag[0], mag[4]]))
    genome_functions = {}
    for ftr in genome["data"]["features"]:
        if "functions" in ftr:
            for func in ftr["functions"]:
                if func not in functions:
                    # display(ftr)
                    ftrid = ftr["id"]+f"_{count}"
                    count += 1
                    result = hashlib.md5(ftr["protein_translation"].encode())
                    md5 = result.hexdigest()
                    result = hashlib.md5(ftr["dna_sequence"].encode())
                    dnamd5 = result.hexdigest()
                    md5_list.append(md5)
                    functions[func] = {"ftrid":ftrid, "probability":1}
                    features[ftrid] = {
                        "aliases": [],
                        "cdss": [
                            ftrid+".CDS"
                        ],
                        "functions":[function],
                        "dna_sequence": ftr["dna_sequence"],
                        "dna_sequence_length": len(ftr["dna_sequence"]),
                        "id": ftrid,
                        "location": [
                            [
                                ftrid+".contig",
                                1,
                                "+",
                                len(ftr["dna_sequence"])
                            ]
                        ],
                        "md5": dnamd5,
                        "ontology_terms": {},
                        "protein_md5": md5,
                        "protein_translation": ftr["protein_translation"],
                        "protein_translation_length": len(ftr["protein_translation"]),
                        "warnings": []
                    }
                    cdsftr = features[ftrid].copy()
                    cdsftr["id"] = ftrid+".CDS"
                    cdsftr["parent_gene"] = ftrid
                    if "aliases" in ftr and len(ftr["aliases"]) >= 1:
                        features[ftrid]["aliases"].append(["gene",ftr["aliases"][0][1]])
                        if ftrid not in feature_aliases:
                            feature_aliases[ftrid] = []
                        feature_aliases[ftrid].append([ftr["aliases"][0][1]])
                elif func not in genome_functions:#Don't want to count same function twice in a genome
                    ftrid = functions[func]["ftrid"]
                    functions[func]["probability"] += 1
                    if "aliases" in ftr and len(ftr["aliases"]) >= 1:
                        features[ftrid]["aliases"].append(["gene",ftr["aliases"][0][1]])
                        if ftrid not in feature_aliases:
                            feature_aliases[ftrid] = []
                        feature_aliases[ftrid].append(ftr["aliases"][0][1])
                genome_functions[func] = True
                functions[func]["probability"] = functions[func]["probability"]/len(mags)
                feature_probabilities[functions[func]["ftrid"]] = functions[func]["probability"]
    md5_list.sort()
    result = hashlib.md5(";".join(md5_list).encode())
    #Writing FASTA
    if not os.path.exists("Assemblies"):
        os.mkdir("Assemblies")
    ofile = open("Assemblies/"+mag[1]+".fasta", "w")
    for func in functions:
        ofile.write(">" + functions[func]["ftrid"] + "\n" +features[functions[func]["ftrid"]]["dna_sequence"] + "\n")
    ofile.close()
    #Saving genome
util.kbdevutil.save("feature_aliases", feature_aliases)
util.kbdevutil.save("functions", functions)
util.kbdevutil.save("feature_probabilities", feature_probabilities)