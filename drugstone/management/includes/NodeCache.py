from collections import defaultdict
import drugstone.models as models


class NodeCache:

    proteins = dict()
    entrez_to_uniprot = defaultdict(lambda: set())
    gene_name_to_uniprot = defaultdict(lambda: set())
    disorders = dict()
    drugs = dict()

    drug_updates = set()
    disorder_updates = set()
    protein_updates = set()

    def init_protein_maps(self):
        print("Generating protein id maps...")
        for protein in self.proteins.values():
            self.entrez_to_uniprot[protein.entrez].add(protein.uniprot_code)
            self.gene_name_to_uniprot[protein.gene].add(protein.uniprot_code)

    def init_proteins(self):
        if len(self.proteins) == 0:
            print("Generating protein maps...")
            for protein in models.Protein.objects.all():
                if protein.id < 1000:
                    protein.delete()
                    continue
                self.proteins[protein.uniprot_code] = protein
        if len(self.proteins) > 0 and (len(self.entrez_to_uniprot) == 0 or len(self.gene_name_to_uniprot) == 0):
            self.init_protein_maps()

    def init_drugs(self):
        if len(self.drugs) == 0:
            print("Generating drug map...")
            for drug in models.Drug.objects.all():
                if drug.id < 1000:
                    drug.delete()
                    continue
                self.drugs[drug.drug_id] = drug

    def init_disorders(self):
        if len(self.disorders) == 0:
            print("Generating disorder map...")
            for disorder in models.Disorder.objects.all():
                if disorder.id < 1000:
                    disorder.delete()
                    continue
                self.disorders[disorder.mondo_id] = disorder

    def is_new_protein(self, protein:models.Protein):
        return protein.uniprot_code in self.protein_updates

    def is_new_drug(self, drug:models.Drug):
        return drug.drug_id in self.drug_updates

    def is_new_disease(self, disease:models.Disorder):
        return disease.mondo_id in self.disorder_updates

    def get_protein_by_uniprot(self,uniprot_id):
        return self.proteins[uniprot_id]

    def get_proteins_by_entrez(self,entrez_id):
        out = list()
        for g in self.entrez_to_uniprot[entrez_id]:
            out.append(self.proteins[g])
        return out

    def get_proteins_by_gene(self, gene_name):
        out = list()
        for g in self.gene_name_to_uniprot[gene_name]:
            out.append(self.proteins[g])
        return out

    def get_drug_by_drugbank(self, drugbank_id):
        return self.drugs[drugbank_id]

    def get_disorder_by_mondo(self, mondo_id):
        return self.disorders[mondo_id]