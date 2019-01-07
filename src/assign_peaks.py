from argparse import ArgumentParser


def main():
    options = get_options()
    gff_data = read_gff(options.ref_gff)
    process_peaks(options, gff_data)


def get_options(arguments=None):
    '''read command line arguments and return the corresponding arguments'''
    parser = ArgumentParser(description='Assign ATACseq peaks to genes')

    parser.add_argument("-i", "--input",
                        dest="peak_centers",
                        help='input peak centers')
    parser.add_argument("-a", "--all_output",
                        dest="all_outputs",
                        help='all features found')
    parser.add_argument("-p", "--priority_output",
                        dest="priority",
                        help='only highest priority features')
    parser.add_argument("-r", "--reference_gff",
                        dest="ref_gff",
                        help='reference general feature file')
    parser.add_argument("-z", "--promoter_size",
                        dest="promoter_size",
                        type=int,
                        help='upstream distance to define promoters')

    return parser.parse_args(arguments)


class Feature():
    '''class for parsing gff lines and returning features'''
    def __init__(self, line):
        tokens = line.strip().split('\t')
        self.scaffold = tokens[0]
        self.type = tokens[2]
        self.start = int(tokens[3])
        self.end = int(tokens[4])
        self.strand = tokens[6]
        self.attributes = tokens[-1]

    def get_id(self):
        self.id = None
        # gene ids are always the first element in the info column
        if self.type == 'gene':
            try:
                self.id = self.attributes.split(';')[0].split('ID=')[1]
            except Exception as e:
                raise(e)

        # transcript ids can be further down
        else:
            for attribute in self.attributes.split(';'):
                token = attribute.split('=')
                if token[0] == 'transcript_id':
                    self.id = token[1]
                    break

    def get_gene(self):
        self.get_id()
        return Gene(self.id,
                    self.start,
                    self.end,
                    self.strand,
                    self.scaffold)


def read_gff(gff_file):
    '''process gff file, building collection of genes keyed by scaffold'''
    with open(gff_file, 'r') as gff_reader:
        last_gene = None
        scaffolds = {}
        for line in gff_reader:
            feature = Feature(line)
            if feature.scaffold not in scaffolds.keys():
                scaffolds[feature.scaffold] = []

            if feature.type == "gene":
                # have to store the previous gene
                if last_gene is not None:
                    scaffolds[last_gene.scaffold].append(last_gene)
                last_gene = feature.get_gene()

            elif feature.type == "exon":
                feature.get_id()
                if feature.id is None:
                    continue
                if last_gene.mrna == feature.id:
                    last_gene.add_feature(feature)

            elif feature.type == "mRNA":
                if last_gene.mrna is None:
                    feature.get_id()
                    if feature.id is None:
                        continue
                    last_gene.mrna = feature.id

            # these all need to be handled separately to find the gene
            elif feature.type in ["five_prime_UTR",
                                  "three_prime_UTR",
                                  "intron"]:

                    # record "last gene",  assuming all introns
                    # are after all genes
                    if last_gene is not None:
                        scaffolds[last_gene.scaffold].append(last_gene)
                        last_gene = None

                    feature.get_id()
                    for gene in scaffolds[feature.scaffold]:
                        if feature.id == gene.mrna:
                            gene.add_feature(feature)
                            break

    return scaffolds


def process_peaks(options, gff_data):
    '''assign each peak to a feature and output the result to all outputs'''
    with open(options.all_outputs, 'w') as all_output, \
            open(options.peak_centers, 'r') as peaks:

        all_output.write("chr\tcoord\tproxtype\tgene\n")
        for line in peaks:
            scaffold, position = line.split()
            for gene, type in process_peak(int(position),
                                           scaffold,
                                           gff_data,
                                           options.promoter_size):
                all_output.write(
                    f"{scaffold}\t{position}\t{type}\t{gene}\n")


def process_peak(position, scaffold, gff_data, promoter_size):
    '''return list of gene, type tuples for the position'''

    if scaffold not in gff_data:
        return [("NA", "intergenic")]

    features = {}
    for gene in gff_data[scaffold]:
        feature = gene.proximity(position, promoter_size)
        if feature is None:
            continue
        features[gene.name] = feature

    # any features added
    if features:
        return [(gene_name, feature)
                for gene_name, feature in features.items()]

    return [("NA", "intergenic")]


class Gene:
    def __init__(self, name, start,
                 end, strand, scaffold):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.scaffold = scaffold
        self.features = {'five_prime_UTR': [],
                         'three_prime_UTR': [],
                         'intron': [],
                         'exon': []}
        self.mrna = None

    def add_feature(self, feature):
        if feature.type == 'five_prime_UTR' or \
                feature.type == 'three_prime_UTR':
            self.features[feature.type].append(
                (feature.start-1,
                 feature.end+1))
        else:
            self.features[feature.type].append(
                (feature.start+1,
                 feature.end-1))

    def coord_in(self, target_coord, target_list):
        for coords in target_list:
            if (target_coord >= coords[0] and target_coord <= coords[1]) or \
                    (target_coord >= coords[1] and target_coord <= coords[0]):
                return True

    def proximity(self, target_coord, promoter_size):
        if target_coord < self.end + 10000 and \
                target_coord > self.start - 10000:

            if self.coord_in(target_coord, self.features['intron']):
                return "intron"
            elif self.coord_in(target_coord, self.features['five_prime_UTR']):
                return "five_utr"
            elif self.coord_in(target_coord, self.features['three_prime_UTR']):
                return "three_utr"
            elif self.coord_in(target_coord, self.features['exon']):
                return "exon"

            elif self.strand == "+":
                if target_coord <= self.start:
                    if target_coord >= self.start - promoter_size:
                        return "promoter"
                    else:
                        return "upstream"
                else:
                    return "downstream"

            elif self.strand == "-":
                if target_coord >= self.end:
                    if target_coord <= self.end + promoter_size:
                        return "promoter"
                    else:
                        return "upstream"
                else:
                    return "downstream"

        else:
            return None


if __name__ == '__main__':
    main()
