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
        dic['target'] = target
        dic['enzyme'] = enzyme
        dic['carrier'] = carrier
        dic['transporter'] = transporter
        return dic


if __name__ == '__main__':
    q = Query()
    print q.query("DB00004")
