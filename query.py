# coding:utf-8
import csv
from collections import defaultdict


class Query:
    def __init__(self):
        self.drugs = self.read("drugs.csv")
        self.partner_protein = self.read("partner_protein.csv")
        self.drug2target = self.read_pair("drug2target.csv")
        self.drug2enzyme = self.read_pair("drug2enzyme.csv")
        self.drug2carrier = self.read_pair("drug2carrier.csv")
        self.drug2transporter = self.read_pair("drug2transporter.csv")
    
    def read(self, csvfile):
        tree = lambda: defaultdict(tree)
        retree = tree()
        ind = 0
        entry = []
        for linelist in csv.reader(open(csvfile)):
            if ind == 0:
                entry = linelist[:]
            else:
                for i in range(1, len(linelist)):
                    retree[linelist[0]][entry[i]] = linelist[i]
            ind += 1
        return retree

    def read_pair(self, csvfile):
        tree = lambda: defaultdict(tree)
        retree = tree()
        ind = 0
        for linelist in csv.reader(open(csvfile)):
            if ind == 0:
                pass
            else:
                drug, temp, act = linelist
                retree[drug][temp] = act
            ind += 1
        return retree

    def query(self, drugid, partner='partner_name'):
        """partner: partner_name, gene_name, uniprot_id"""
        drugname = self.drugs[drugid]['drugname']
        drugtype = self.drugs[drugid]['drug_type']
        drugATC = self.drugs[drugid]['ATC_codes']
        drugA = self.drugs[drugid]['approved']
        target = [self.partner_protein[i][partner]
                  for i, j in self.drug2target[drugid].items()]
        enzyme = [self.partner_protein[i][partner]
                  for i, j in self.drug2enzyme[drugid].items()]
        carrier = [self.partner_protein[i][partner]
                   for i, j in self.drug2carrier[drugid].items()]
        transporter = [self.partner_protein[i][partner]
                       for i, j in self.drug2transporter[drugid].items()]
        dic = defaultdict(dict)
        dic['drugname'] = drugname
        dic['drugtype'] = drugtype
        dic['drugATC'] = drugATC
        dic['target'] = target
        dic['enzyme'] = enzyme
        dic['carrier'] = carrier
        dic['transporter'] = transporter
        dic['approved'] = drugA
        return dic

    def dotfile(self, drugids, dotfile):
        graphformat = open("template.dot").read()
        nodes = []
        edges = []
        for drugid in drugids:
            dic = self.query(drugid)
            nodes.append((dic['drugname'], dic['drugtype']))
            for terms in ['target', 'enzyme', 'carrier', 'transporter']:
                for term in dic[terms]:
                    nodes.append((term, terms))
                    edges.append((dic['drugname'], term))
        res = ([], [], [], [], [], [])
        for node in nodes:
            try:
                ind = ['target', 'enzyme',
                       'carrier', 'transporter'].index(node[1])
                res[ind+1].append('"' + node[0] + '"')
            except ValueError:
                res[0].append('"' + node[0] + '"')
        for edge in edges:
            res[5].append('"{}" -- "{}"'.format(edge[0], edge[1]))

        content = map(lambda t: "\n    " + ";\n    ".join(t) + ";\n", res)
        with open(dotfile, 'w') as fw:
            fw.write(graphformat % tuple(content))

    def tofile(self, drugids, filename):
        csvfile = open(filename, 'wb')
        writer = csv.writer(csvfile)
        writer.writerow(["drugID", "drugname", "drugtype", "targets",
                         "enzymes", "carriers", "transporters"])
        for drugid in drugids:
            dic = self.query(drugid)
            if dic['approved'] == '1':
                writer.writerow([drugid, dic['drugname'], dic['drugtype'],
                                 dic['target'], dic['enzyme'], dic['carrier'],
                                 dic['transporter']])
        csvfile.close()


if __name__ == '__main__':
    q = Query()
    q.tofile(q.drugs.keys(), 'total.csv')
    # q.dotfile(['DB00004', 'DB00005'], dotfile="graph.dot")
