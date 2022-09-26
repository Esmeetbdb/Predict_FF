import numpy as np
import pickle

def coverage(contig_dictionary,bam_file,bin_size,min_qual,max_len):
    import pysam
    import count_reads

    samfile = pysam.AlignmentFile(bam_file, "r")
    insert_size_per_region = {}
    coverage_per_region = {}

    for contig in contig_dictionary:
        coverage_per_region[contig] = []
        for i in range(0,contig_dictionary[contig],bin_size):
            coverage_per_region[contig].append(0)
            for read in samfile.fetch(contig, i, (i+bin_size)):
                read_length = read.template_length
                if read_length >= 50 and count_reads.read_ok(read, min_qual) and read_length <= max_len:
                    coverage_per_region[contig][-1] += 1

    autosomal=[]
    for contig in coverage_per_region:	
        if "X" in contig or "Y" in contig:
           continue
        autosomal+=coverage_per_region[contig]
    med_autosomal=np.median(autosomal)

    if "chrY" in coverage_per_region:
        Y="chrY"
        X="chrX"
    else:
        Y="Y"
        X="X"

    print("#chromosome\tstart\tend\fragment_fraction")
    for i in range(0,len(coverage_per_region[Y])):
         print("{}\t{}\t{}\t{}".format(Y,(i-1)*bin_size+1,i*bin_size,coverage_per_region[Y][i]/med_autosomal ) )

    for i in range(0,len(coverage_per_region[X])):
         print("{}\t{}\t{}\t{}".format(X,(i-1)*bin_size+1,i*bin_size,coverage_per_region[X][i]/med_autosomal ) )

    return()

def summarise(samplesheet):
	header=["Individual","FFY"]
	header_printed=False
	
	for f in open(samplesheet):
		sample=f.split()[0].split("/")[-1].split(".")[0]
		ff=f.split()[1]		
		out=[sample,ff]

		for line in open( f.split()[0] ):
			if line[0] == "#":
				continue
			content=line.strip().split()
			if not header_printed:
				header.append(content[0])
			out.append(content[-1])

		if not header_printed:
			print(",".join(header))
			header_printed=True

		print(",".join(out))

def retrieve_ffy_model(v,model_path):
    import pandas
    v=pandas.DataFrame(data=v, index=[0])
    
    with open(model_path, 'rb') as file:
        data = pickle.load(file)

    pca = data["pca"]
    model=data["model"]
    scaler = data["scaler"]
    predictors = data["predictors"]
    v=v[predictors].values
    
    x_scaled = scaler.transform(v.reshape(1, -1))
    x_scaled = pca.transform(x_scaled)
    res = model.predict(x_scaled)

    return(res)
