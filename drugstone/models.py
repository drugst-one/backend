from django.core.exceptions import ValidationError
from django.db import models


# Main biological and medical entities


class PPIDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'
    
    class Meta:
        unique_together = ("name", "version", "licenced")


class PDIDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'

    class Meta:
        unique_together = ("name", "version", "licenced")


class PDisDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'

    class Meta:
        unique_together = ("name", "version", "licenced")


class DrDiDataset(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=False)
    link = models.CharField(max_length=128, default="", unique=False)
    version = models.CharField(max_length=128, default="", unique=False)
    licenced = models.BooleanField(default=True)

    def __str__(self):
        return f'{self.name}-{self.version}_{"licenced" if self.licenced else "unlicenced"}'

    class Meta:
        unique_together = ("name", "version", "licenced")


class EnsemblGene(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=15)  # starts with ENSG...
    protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, related_name="ensg"
    )
    
class CellularComponent(models.Model):
    id = models.AutoField(primary_key=True)
    go_code = models.CharField(max_length=10)
    display_name = models.CharField(max_length=255, default="")
    layer = models.CharField(max_length=255, default="other")
    
    class Meta:
        unique_together = ("go_code", "display_name")
    
    def __str__(self):
        return self.display_name

    def __eq__(self, other):
        return (
            self.go_code == other.go_code
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.go_code))

    def update(self, other):
        self.go_code = other.go_code
        self.display_name = other.display_name
        self.description = other.description
        self.layer = other.layer


class Protein(models.Model):
    # According to https://www.uniprot.org/help/accession_numbers UniProt accession codes
    # are either 6 or 10 characters long
    id = models.AutoField(primary_key=True)
    uniprot_code = models.CharField(max_length=10)
    gene = models.CharField(max_length=127, default="")  # symbol
    protein_name = models.CharField(max_length=255, default="")
    entrez = models.CharField(max_length=15, default="")
    drugs = models.ManyToManyField(
        "Drug", through="ProteinDrugInteraction", related_name="interacting_drugs"
    )
    cellular_components = models.ManyToManyField(
        "CellularComponent", through="ActiveIn", related_name="active_in"
    )
    tissue_expression = models.ManyToManyField(
        "Tissue", through="ExpressionLevel", related_name="interacting_drugs"
    )

    class Meta:
        unique_together = ("uniprot_code", "gene", "entrez")

    def __str__(self):
        return self.gene

    def __eq__(self, other):
        return (
            self.uniprot_code == other.uniprot_code
            and self.gene == other.gene
            and self.protein_name == other.protein_name
            and self.entrez == other.entrez
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.uniprot_code, self.gene, self.entrez))

    def update(self, other):
        self.uniprot_code = other.uniprot_code
        self.gene = other.gene
        self.protein_name = other.protein_name
        self.entrez = other.entrez


class ExpressionLevel(models.Model):
    id = models.AutoField(primary_key=True)
    tissue = models.ForeignKey("Tissue", on_delete=models.CASCADE)
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE)
    expression_level = models.FloatField()

    class Meta:
        unique_together = ("tissue", "protein")

    def __hash__(self):
        return hash(f"{self.tissue_id}_{self.protein_id}")

class ActiveIn(models.Model):
    id = models.AutoField(primary_key=True)
    cellularComponent = models.ForeignKey("CellularComponent", on_delete=models.CASCADE)
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE)
    class Meta:
        unique_together = ("cellularComponent", "protein")
    def __hash__(self):
        return hash((self.cellularComponent_id, self.protein_id))

class Tissue(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, default="", unique=True)

    def __str__(self):
        return self.name


class Disorder(models.Model):
    id = models.AutoField(primary_key=True)
    mondo_id = models.CharField(max_length=7)
    label = models.CharField(max_length=256, default="")  # symbol
    icd10 = models.CharField(max_length=512, default="")
    proteins = models.ManyToManyField(
        "Protein",
        through="ProteinDisorderAssociation",
        related_name="associated_proteins",
    )

    class Meta:
        unique_together = ("mondo_id", "label", "icd10")

    def __str__(self):
        return self.label

    def __eq__(self, other):
        return (
            self.mondo_id == other.mondo_id
            and self.label == other.label
            and self.icd10 == other.icd10
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.mondo_id, self.label, self.icd10))

    def update(self, other):
        self.mondo_id = other.mondo_id
        self.label = other.label
        self.icd10 = other.icd10


class Drug(models.Model):
    id = models.AutoField(primary_key=True)
    drug_id = models.CharField(max_length=10, unique=True)
    name = models.CharField(max_length=256, default="")
    status = models.CharField(max_length=128, default="")
    # in_trial = models.BooleanField(default=False)
    # in_literature = models.BooleanField(default=False)
    links = models.CharField(max_length=16 * 1024, default="")

    def __str__(self):
        return self.drug_id

    def __eq__(self, other):
        return (
            self.drug_id == other.drug_id
            and self.name == other.name
            and self.status == other.status
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.drug_id, self.name, self.status))

    def update(self, other):
        self.drug_id = other.drug_id
        self.name = other.name
        self.status = other.status
        self.links = other.links


class ProteinDisorderAssociation(models.Model):
    id = models.BigAutoField(primary_key=True)
    pdis_dataset = models.ForeignKey(
        "PDisDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="pdis_dataset_relation",
    )
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE)
    disorder = models.ForeignKey("Disorder", on_delete=models.CASCADE)
    score = models.FloatField()

    class Meta:
        unique_together = ("pdis_dataset", "protein", "disorder")

    def __str__(self):
        return f"{self.pdis_dataset}-{self.protein}-{self.disorder}"

    def __eq__(self, other):
        return (
            self.pdis_dataset_id == other.pdis_dataset_id
            and self.protein_id == other.protein_id
            and self.disorder_id == other.disorder_id
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.pdis_dataset_id, self.protein_id, self.disorder_id))


class DrugDisorderIndication(models.Model):
    id = models.AutoField(primary_key=True)
    drdi_dataset = models.ForeignKey(
        "DrDiDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="drdi_dataset_relation",
    )
    drug = models.ForeignKey("Drug", on_delete=models.CASCADE)
    disorder = models.ForeignKey("Disorder", on_delete=models.CASCADE)

    class Meta:
        unique_together = ("drdi_dataset", "drug", "disorder")

    def __str__(self):
        return f"{self.drdi_dataset}-{self.drug}-{self.disorder}"

    def __eq__(self, other):
        return (
            self.drdi_dataset_id == other.drdi_dataset_id
            and self.drug_id == other.drug_id
            and self.disorder_id == other.disorder_id
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.drdi_dataset_id, self.drug_id, self.disorder_id))


class ProteinProteinInteraction(models.Model):
    id = models.BigAutoField(primary_key=True)
    ppi_dataset = models.ForeignKey(
        "PPIDataset",
        null=True,
        on_delete=models.CASCADE,
        related_name="ppi_dataset_relation",
    )
    from_protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, related_name="interacting_proteins_out"
    )
    to_protein = models.ForeignKey(
        "Protein", on_delete=models.CASCADE, related_name="interacting_proteins_in"
    )

    def validate_unique(self, exclude=None):
        p1p2_q = ProteinProteinInteraction.objects.filter(
            from_protein=self.from_protein,
            to_protein=self.to_protein,
            ppi_dataset=self.ppi_dataset,
        )
        p2p1_q = ProteinProteinInteraction.objects.filter(
            from_protein=self.to_protein,
            to_protein=self.from_protein,
            ppi_dataset=self.ppi_dataset,
        )

        if p1p2_q.exists() or p2p1_q.exists():
            raise ValidationError("Protein-Protein interaction must be unique!")

    def save(self, *args, **kwargs):
        self.validate_unique()
        super(ProteinProteinInteraction, self).save(*args, **kwargs)

    def __str__(self):
        return f"{self.ppi_dataset}-{self.from_protein}-{self.to_protein}"

    def __eq__(self, other):
        return (
            self.ppi_dataset_id == other.ppi_dataset_id
            and self.from_protein_id == other.from_protein_id
            and self.to_protein_id == other.to_protein_id
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.ppi_dataset_id, self.from_protein_id, self.to_protein_id))
    

class ProteinDrugInteraction(models.Model):
    id = models.BigAutoField(primary_key=True)
    pdi_dataset = models.ForeignKey(
        PDIDataset,
        null=True,
        on_delete=models.CASCADE,
        related_name="pdi_dataset_relation",
    )
    protein = models.ForeignKey("Protein", on_delete=models.CASCADE)
    drug = models.ForeignKey("Drug", on_delete=models.CASCADE)
    actions = models.CharField(max_length=255, default="[]")

    class Meta:
        unique_together = ("pdi_dataset", "protein", "drug")

    def __str__(self):
        return f"{self.pdi_dataset}-{self.protein}-{self.drug}"

    def __eq__(self, other):
        return (
            self.pdi_dataset_id == other.pdi_dataset_id
            and self.protein_id == other.protein_id
            and self.drug_id == other.drug_id
            and self.actions == other.actions
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.pdi_dataset_id, self.protein_id, self.drug_id, self.actions))
    
    # def bulk_create(self, objs, existing: set, **kwargs):
    #     history_objs = create_history_objs(objs, existing)
    #     # store changes in history
    #     models.ProteinDrugInteractionHistory.bulk_create(history_objs, **kwargs)
    #     return super(ProteinDrugInteraction, self).bulk_create(objs, **kwargs)
    
    # def historic(self, objs, historic_dataset: PPIDataset):
    #     # load historic version of objs
    #     historic_changes = models.ProteinDrugInteractionHistory.filter(pdi_dataset__name=historic_dataset.name)
    #     print('historic_changes', historic_changes)
    #     # convert list of objs to dict for faster lookup, key is composed out of unique identifier for element
    #     objs = left_search(objs, historic_changes, historic_dataset)
    #     return objs

# def left_search(objs, historic_changes, historic_dataset, key_attributes):
#     objs_dict = {hash((getattr(change, attr) for attr in key_attributes)): obj for obj in objs}
#     left_search_done = set()
#     # changes are sorted by version, hence in correct order
#     for change in historic_changes:
#         key = hash((getattr(change, attr) for attr in key_attributes))
#         if key in left_search_done:
#             # left-search is done for the element affected by this change
#             continue
        
#         if change.new and change.pdi_dataset.id > historic_dataset.id:
#             # item was added in this version, hence it did not exist in previous version
#             del objs_dict[key]
#         elif change.pdi_dataset.id <= historic_dataset.id:
#             # update and return, left-search terminated, state of target dataset reached
#             objs_dict[key] = change
#             left_search_done.add(key)
#         else:
#             # change.new --> False && change.pdi_dataset.id <= historic_dataset.id
#             # update but continue left search
#             objs_dict[key] = change
#     return objs_dict.values()


class Task(models.Model):
    token = models.CharField(max_length=32, unique=True, primary_key=True)
    created_at = models.DateTimeField(auto_now_add=True)
    target = models.CharField(
        max_length=32, choices=[("drug", "Drug"), ("drug-target", "Drug Target")]
    )

    algorithm = models.CharField(max_length=128)
    parameters = models.TextField()

    progress = models.FloatField(default=0.0)  # Progress as fraction (0.0 - 1.0)
    started_at = models.DateTimeField(null=True)
    finished_at = models.DateTimeField(null=True)
    worker_id = models.CharField(max_length=128, null=True)
    job_id = models.CharField(max_length=128, null=True)
    done = models.BooleanField(default=False)
    failed = models.BooleanField(default=False)
    status = models.CharField(max_length=255, null=True)

    result = models.TextField(null=True)


class Network(models.Model):
    id = models.CharField(primary_key=True, max_length=32, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    nodes = models.TextField(null=True, default="")
    edges = models.TextField(null=True, default="")
    config = models.TextField(null=True, default="")
    groups = models.TextField(null=True, default="")


#### save history for versioning ####

# def create_history_objs(objs, existing):
#     history_objs = objs.copy()
#     for hobj in history_objs:
#         hobj.new = hobj not in existing
#     return history_objs

# class EnsemblGeneHistory(EnsemblGene):
#     version = models.CharField(max_length=32)

 
# class ProteinHistory(Protein):
#     version = models.CharField(max_length=32)

    
# class DisorderHistory(Disorder):
#     version = models.CharField(max_length=32)


# class DrugHistory(Drug):
#     version = models.CharField(max_length=32)

    
# class ProteinProteinInteractionHistory(ProteinProteinInteraction):
#     class Meta:
#         ordering = ('ppi_dataset.name', '-ppi_dataset.id')
    
    
# class ProteinDisorderAssociationHistory(ProteinDisorderAssociation):    
#     class Meta:
#         ordering = ('pdis_dataset.name', '-pdis_dataset.id')

    
# class ProteinDrugInteractionHistory(ProteinDrugInteraction):
#     new = models.BooleanField(default=True)
#     class Meta:
#         ordering = ('pdi_dataset.name', '-pdi_dataset.id')

