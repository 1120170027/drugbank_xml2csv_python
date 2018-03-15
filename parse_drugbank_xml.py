# coding:utf-8
import csv
from collections import defaultdict
from StringIO import StringIO

from lxml import etree
import sys
reload(sys)
sys.setdefaultencoding("utf-8")

# Parse XML
f = open('drugbank.xml', 'r')
data = f.read()
f.close()

tree = etree.parse(StringIO(data))
context = etree.iterparse(StringIO(data))

root = tree.getroot()
# print len(root), 'drugs'

#######################################################################
# Iterate over drugs
drug2attrib = defaultdict(dict)
# drugbank_id -> {'drugname', 'drug_type', 'groups',
# 'targets/enzymes/transporters': [_id, _actions]}

template2attrib = {
    'targets': defaultdict(dict),
    'enzymes': defaultdict(dict),
    'transporters': defaultdict(dict),
    'carriers': defaultdict(dict)
}

sin = {
    "targets": "target",
    "enzymes": "enzyme",
    "transporters": "transporter",
    "carriers": "carrier"
}
# drugbank_target_id -> {'gene', 'name', 'organism', 'taxonomy_id',
# 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}

tag_prefix = '{http://www.drugbank.ca}'

fw = open("description.csv", 'w')
writer = csv.writer(fw)

for child in root:
    for s in child.findall(tag_prefix + 'drugbank-id'):
        if 'primary' in s.attrib:
            drugbank_id = s.text

    drugname = child.findall(tag_prefix + 'name')[0].text
    drug2attrib[drugbank_id]['drugname'] = drugname

    drug_type = child.attrib['type']
    drug2attrib[drugbank_id]['drug_type'] = drug_type

    groups = [
        s.text
        for s in child.find(tag_prefix + 'groups').findall(
            tag_prefix + 'group')
    ]
    drug2attrib[drugbank_id]['groups'] = groups

    if child.findall(tag_prefix + 'atc-codes')[0].text is None:
        drug2attrib[drugbank_id]['ATC_codes'] = 0
    else:
        drug2attrib[drugbank_id]['ATC_codes'] = 1

    writer.writerow([drugbank_id, drugname, child.find(tag_prefix + 'description').text])

    # Get targets, enzymes, transporters, carriers
    for template in ['targets', 'enzymes', 'transporters', 'carriers']:
        drug2attrib[drugbank_id][template] = []
        for res in child.find(tag_prefix + template).findall(
                tag_prefix + sin[template]):
            if res.find(tag_prefix + 'polypeptide') is None:
                continue
            template_id = res.find(tag_prefix + 'id').text
            template_gene = res.find(tag_prefix + 'polypeptide').find(
                tag_prefix + 'gene-name').text
            template_name = res.find(tag_prefix + 'name').text
            template_organism = res.find(tag_prefix + 'organism').text
            template_taxonomy_id = res.find(tag_prefix + 'polypeptide').find(
                tag_prefix + 'organism').attrib['ncbi-taxonomy-id']

            if template_organism is None and template_taxonomy_id == '9606':
                template_organism = 'Human'
            if template_organism == 'Human' and template_taxonomy_id == '':
                template_taxonomy_id = '9606'
            if template_organism == 'Homo sapiens':
                template_organism = 'Human'
            if template_gene is None or template_organism is None \
                    or template_taxonomy_id is None:
                continue

            template_external_ids = res.find(tag_prefix + 'polypeptide').find(
                tag_prefix + 'external-identifiers').findall(
                    tag_prefix + 'external-identifier')
            template_uniprot_id = ''
            template_genbank_gene = ''
            template_genbank_protein = ''
            template_hgnc_id = ''

            for external_id in template_external_ids:
                if external_id.find(
                        tag_prefix + 'resource').text == 'UniProtKB':
                    template_uniprot_id = external_id.find(
                        tag_prefix + 'identifier').text
                elif external_id.find(tag_prefix + 'resource'
                                      ).text == 'GenBank Gene Database':
                    template_genbank_gene = external_id.find(
                        tag_prefix + 'identifier').text
                elif external_id.find(tag_prefix + 'resource'
                                      ).text == 'GenBank Protein Database':
                    template_genbank_protein = external_id.find(
                        tag_prefix + 'identifier').text
                elif external_id.find(
                        tag_prefix + 'resource'
                ).text == 'HUGO Gene Nomenclature Committee (HGNC)':
                    template_hgnc_id = external_id.find(
                        tag_prefix + 'identifier').text

            template_actions = [
                s.text.lower()
                for s in res.find(tag_prefix + 'actions').findall(
                    tag_prefix + 'action')
            ]

            drug2attrib[drugbank_id][template].append((template_id,
                                                       template_actions))

            # {'gene', 'name', 'organism', 'taxonomy_id', 'uniprot_id',
            # 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}
            if template_id not in template2attrib[template]:
                template2attrib[template][template_id]['gene'] = template_gene
                template2attrib[template][template_id]['name'] = template_name
                template2attrib[template][template_id][
                    'organism'] = template_organism
                template2attrib[template][template_id][
                    'taxonomy_id'] = template_taxonomy_id

                template2attrib[template][template_id][
                    'uniprot_id'] = template_uniprot_id
                template2attrib[template][template_id][
                    'genbank_gene_id'] = template_genbank_gene
                template2attrib[template][template_id][
                    'genbank_protein_id'] = template_genbank_protein
                template2attrib[template][template_id][
                    'hgnc_id'] = template_hgnc_id

fw.close()
print '\n'

#######################################################################
# List of drugs to save (as long as num_targets + num_enzymes +
# num_transporters != 0)
drugs = []
for drugbank_id in sorted(drug2attrib.keys()):
    if len(drug2attrib[drugbank_id]['targets']) == 0 and \
            len(drug2attrib[drugbank_id]['enzymes']) == 0 and \
            len(drug2attrib[drugbank_id]['transporters']) == 0 and \
            len(drug2attrib[drugbank_id]['carriers']) == 0:
        continue
    else:
        drugs.append(drugbank_id)

print len(drug2attrib), "drugs parsed from XML"
print len(
    drugs), "drugs with at least 1 target/ enzyme/ transporter / carriers"

#######################################################################
# Save drug attributes to CSV {'drugname', 'drug_type', 'groups',
# 'targets/enzymes/transporters': [_id, _actions]}
outf = open('drugs.csv', 'w')
writer = csv.writer(outf)
writer.writerow([
    'drugbank_id', 'drugname', 'drug_type', 'ATC_codes', 'approved',
    'experimental', 'illicit', 'investigational', 'nutraceutical', 'withdrawn'
])

for drugbank_id in drugs:
    drugname = drug2attrib[drugbank_id]['drugname']
    if isinstance(drugname, unicode):
        if u'\u03b2' in drugname:
            drugname = drugname.replace(u'\u03b2', 'beta')
        if u'\u03b1' in drugname:
            drugname = drugname.replace(u'\u03b1', 'alpha')
        drugname = drugname.encode("utf-8")
    drug_type = drug2attrib[drugbank_id]['drug_type']
    ATC_codes = 1 if drug2attrib[drugbank_id]['ATC_codes'] == 1 else 0
    groups = [
        1 if group in drug2attrib[drugbank_id]['groups'] else 0
        for group in [
            'approved', 'experimental', 'illicit', 'investigational',
            'nutraceutical', 'withdrawn'
        ]
    ]

    writer.writerow([drugbank_id, drugname, drug_type, ATC_codes] + groups)

outf.close()

#######################################################################
# Save drug-target, -enzyme, -transporter pairs to CSV

# target [('antagonist', 1374), ('agonist', 857), ('inhibitor', 1818)]
# enzyme [('substrate', 2402), ('inducer', 407), ('inhibitor', 1350)]
# transporter [('substrate', 790), ('inducer', 100), ('inhibitor', 1075)]

# Targets
outf = open('drug2target.csv', 'w')
outfh = open('drug2target_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

writer.writerow(['drugbank_id', 'partner_id', 'actions'])
writerh.writerow(['drugbank_id', 'partner_id', 'actions'])

for drugbank_id in drugs:
    for (target_id, target_actions) in drug2attrib[drugbank_id]['targets']:
        target_actions = "#".join(target_actions)
        writer.writerow([drugbank_id, target_id, target_actions])
        if template2attrib['targets'][target_id]['organism'] == 'Human' and \
                template2attrib['targets'][target_id]['taxonomy_id'] == '9606':
            writerh.writerow([drugbank_id, target_id, target_actions])

outf.close()
outfh.close()

# Enzymes
outf = open('drug2enzyme.csv', 'w')
outfh = open('drug2enzyme_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

writer.writerow(['drugbank_id', 'partner_id', 'actions'])
writerh.writerow(['drugbank_id', 'partner_id', 'actions'])

for drugbank_id in drugs:
    for (enzyme_id, enzyme_actions) in drug2attrib[drugbank_id]['enzymes']:
        enzyme_actions = "#".join(enzyme_actions)
        writer.writerow([drugbank_id, enzyme_id, enzyme_actions])
        if template2attrib['enzymes'][enzyme_id]['organism'] == 'Human' and \
                template2attrib['enzymes'][enzyme_id]['taxonomy_id'] == '9606':
            writerh.writerow([drugbank_id, enzyme_id, enzyme_actions])

outf.close()
outfh.close()

# Transporters
outf = open('drug2transporter.csv', 'w')
outfh = open('drug2transporter_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

writer.writerow(['drugbank_id', 'partner_id', 'actions'])
writerh.writerow(['drugbank_id', 'partner_id', 'actions'])

for drugbank_id in drugs:
    for (transporter_id,
         transporter_actions) in drug2attrib[drugbank_id]['transporters']:
        transporter_actions = "#".join(transporter_actions)
        writer.writerow([drugbank_id, transporter_id, transporter_actions])
        if template2attrib['transporters'][transporter_id]['organism'] ==\
                'Human' and\
                template2attrib['transporters'][transporter_id]['taxonomy_id']\
                == '9606':
            writerh.writerow(
                [drugbank_id, transporter_id, transporter_actions])

outf.close()
outfh.close()

# Carriers
outf = open('drug2carrier.csv', 'w')
outfh = open('drug2carrier_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

writer.writerow(['drugbank_id', 'partner_id', 'actions'])
writerh.writerow(['drugbank_id', 'partner_id', 'actions'])

for drugbank_id in drugs:
    for (carrier_id, carrier_actions) in drug2attrib[drugbank_id]['carriers']:
        carrier_actions = "#".join(carrier_actions)
        writer.writerow([drugbank_id, carrier_id, carrier_actions])

        if template2attrib['carriers'][carrier_id]['organism'] == 'Human' \
                and template2attrib['carriers'][carrier_id]['taxonomy_id'] \
                == '9606':
            writerh.writerow([drugbank_id, carrier_id, carrier_actions])

outf.close()
outfh.close()

#######################################################################
# Save all targets, enzymes, transporters, carriers to CSV

# drugbank_target_id -> #{'gene', 'name', 'organism', 'taxonomy_id',
# 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}

outf = open('partner_protein.csv', 'w')
outfh = open('partner_protein_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

writer.writerow([
    'partner_id', 'partner_name', 'gene_name', 'uniprot_id', 'genbank_gene_id',
    'genbank_protein_id', 'hgnc_id', 'organism', 'taxonomy_id'
])
writerh.writerow([
    'partner_id', 'partner_name', 'gene_name', 'uniprot_id', 'genbank_gene_id',
    'genbank_protein_id', 'hgnc_id', 'organism', 'taxonomy_id'
])

partners_written = set()

for template in ['targets', 'enzymes', 'transporters', 'carriers']:
    for partner_id in sorted(template2attrib[template].keys()):
        if partner_id in partners_written:
            # print partner_id, template2attrib[template][partner_id]['gene'],
            # 'already recorded'
            continue
        partner_name = template2attrib[template][partner_id]['name']
        gene_name = template2attrib[template][partner_id]['gene']
        organism = template2attrib[template][partner_id]['organism']
        taxonomy_id = template2attrib[template][partner_id]['taxonomy_id']

        uniprot_id = template2attrib[template][partner_id]['uniprot_id']
        genbank_gene_id = template2attrib[template][partner_id][
            'genbank_gene_id']
        genbank_protein_id = template2attrib[template][partner_id][
            'genbank_protein_id']
        hgnc_id = template2attrib[template][partner_id]['hgnc_id']

        partners_written.add(partner_id)

        writer.writerow([
            partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id,
            genbank_protein_id, hgnc_id, organism, taxonomy_id
        ])

        if taxonomy_id == '9606' and organism == 'Human':
            writerh.writerow([
                partner_id, partner_name, gene_name, uniprot_id,
                genbank_gene_id, genbank_protein_id, hgnc_id, organism,
                taxonomy_id
            ])

        if taxonomy_id == '9606' and organism.lower() != 'human':
            print partner_id, template2attrib[template][partner_id][
                'gene'], organism, taxonomy_id, 'organism mismatch'

outf.close()
outfh.close()
