#!/usr/bin/python

import os
import sys
from py2cytoscape.cyrest import cyclient
import pybiomart as pbm
import requests
from urllib.parse import urlencode, urljoin

variant_file = sys.argv[1]
network_name = str(sys.argv[2])
#network_name = 'BIGCAT task'

# Step 1
# connecting to biomart hsapiens variant database
dataset = pbm.Dataset(name='hsapiens_snp',host='http://www.ensembl.org')
variants = []

# reading the variants from file
with open(variant_file, 'r') as f:
    for line in f:
        variants.append(line.strip())

# retrieving the gene ids by variant ids using biomart attributes
genes_by_var = dataset.query(attributes=['refsnp_id', 'ensembl_gene_stable_id'],
                             filters={'snp_filter': variants})
#print (dataset.attributes)

# Step 2
# create cytoscape client (make sure app is running)
cy = cyclient()
# create a new empty network
cy.network.create_empty(name=network_name)

# add variants to the network as nodes
for snp in variants:
    cy.network.add_node(network=network_name, name=snp)

# getting keys (names) of the two columns
k_snp, k_gen = genes_by_var.keys()
# add genes given by biomart taking unique gene names
for gene in set(genes_by_var[k_gen]):
    cy.network.add_node(network=network_name, name=gene)

# add edges to the network using the dataframe given by biomart query
for index, row in genes_by_var.iterrows():
    cy.network.add_edge(network=network_name, sourceName=row[k_snp], targetName=row[k_gen])

#cy.network.show()

# Step 3
# setting the url and parameters of the cytargetlinker app
ctl_url = 'http://localhost:1234/v1/commands/cytargetlinker/'
# parameters for the extend request using as wikipathways a file in current directory
extend_parameters = {'idAttribute': 'shared name',
                     'linkSetDirectory': os.getcwd(), 'network': network_name}

# sending the request to extend network
r = requests.get(urljoin(ctl_url, "extend?{}").format(urlencode(extend_parameters)))

#print(urljoin(ctl_url, "extend?{}").format(urlencode(extend_parameters)))
#print(r.text)

# updating network name since after extend request its changed by ctl
network_name = "CTL_" + network_name

# apply ctl layout
r = requests.get(urljoin(ctl_url, "applyLayout?{}").format(urlencode({'network': network_name})))
#print("apply layout\n", r.text)

# apply ctl visual style
requests.get(urljoin(ctl_url, "applyVisualstyle?{}").format(urlencode({'network': network_name})))

# save as image
cy.session.save_as(os.path.join(os.getcwd(), "bigcat_task"))
outputFileName = os.path.join(os.getcwd(), "bigcat_task_img.svg")
cy.view.export(options='SVG', outputFile=outputFileName, Resolution='600')
