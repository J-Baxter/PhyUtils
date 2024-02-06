# Edits a beauti XML file to remove clock and tree models to fit trait analyses
# to empirical trees.
#
# Arguments: 1) beauti XML that includes all your trait analysis
#            2) file path to posterior distribution of phylogenetic trees
#            (as obtained from BEAST)
#
# This script works reasonably well on typical BEAST parameterisations and model
# selection. If you have included advanced models in your Beauti XML, you may
# need to supplement 'elements.py'

# Copyright (c) 2024 James Baxter

# import modules
import xml.etree.ElementTree as ET
import argparse
import re
#import keys

keys = [['', 'id', 'alignment'],
        ['', 'id', 'treeModel'],
        ['', 'id', 'treeLength'],
        ['', 'id', 'CP1+2.patterns'],
        ['', 'id', 'CP3.patterns'],
        ['', 'id', 'initialDemo'],
        ['', 'id', 'age(root)'],
        ['', 'id', 'skygrid'],
        ['', 'id', 'branchRates'],
        ['', 'id', 'meanRate'],
        ['', 'id', 'coefficientOfVariation'],
        ['', 'id', 'covariance'],
        ['', 'id', 'CP1+2.hky'],
        ['', 'id', 'CP3.hky'],
        ['', 'id', 'CP1+2.siteModel'],
        ['', 'id', 'CP3.siteModel'],
        ['', 'id', 'allMus'],
        ['', 'id', 'treeLikelihood'],
        ['', 'id', "startingTree"],
        ['', 'id', "allNus"],
        ['parameter', 'idref', 'CP1+2.kappa'],
        ['parameter', 'idref', 'CP3.kappa'],
        ['parameter', 'idref', 'frequencies'],
        ['parameter', 'idref', 'CP1+2.alpha'],
        ['parameter', 'idref', 'CP3.alpha'],
        ['parameter', 'idref', 'ucld.mean'],
        ['parameter', 'idref', 'ucld.stdev'],
        ['parameter', 'idref', 'treeModel.allInternalNodeHeights'],
        ['parameter', 'idref', 'ucld.mean'],
        ['parameter', 'idref', 'branchRates.categories'],
        ['parameter', 'idref', 'allMus'],
        ['parameter', 'idref', 'treeModel.rootHeight'],
        ['parameter', 'idref', 'treeModel.internalNodeHeights'],
        ['parameter', 'idref', 'skygrid.precision'],
        ['parameter', 'idref', 'ucld.mean'],
        ['parameter', 'idref', 'ucld.stdev'],
        ['parameter', 'idref', 'skygrid.logPopSize'],
        ['parameter', 'idref', 'skygrid.cutOff'],
        ['parameter', 'idref', 'skygrid'],
        ['parameter', 'idref', "age(root)"],
        ['parameter', 'idref', "allNus"],
        ['parameter', 'idref', "constant.popSize"],
        ['tmrcaStatistic', 'idref', "age(root)"],
        ['discretizedBranchRates', 'idref', "branchRates"],
        ['treeModel', 'idref', 'treeModel'],
        ['siteModel', 'idref', "CP1+2.siteModel"],
        ['siteModel', 'idref', "CP3.siteModel"],
        ['joint', 'idref', 'joint'],
        ['upDownOperator', '', ''],
        ['subtreeSlide', '', ''],
        ['narrowExchange', '', ''],
        ['wideExchange', '', ''],
        ['wilsonBalding', '', ''],
        ['gmrfGridBlockUpdateOperator', '', ''],
        ['trancated', '', ''],
        ['gmrfSkyGridLikelihood', '', ''],
        ['treeDataLikelihood', '', ''],
        ['discretizedBranchRates', '', ''],
        ['constantSize', '', ''],
        ['coalescentLikelihood', '', ''],
        ['compoundParameter', '', ''],

        ]


def parse_args():
    parser = argparse.ArgumentParser(description="a script to do stuff")
    parser.add_argument("xml_path")
    parser.add_argument("tree_path")
    args = parser.parse_args()
    return args


def parse_query(key):
    q = ['', '', '']

    if key[1] and key[2] != '':
        q[1] = f"[@{key[1]}="

        if all(x != "" for x in key):
            q[2] = f"'{key[2]}']..."
        else:
            q[2] = f"'{key[2]}']"
    else:
        q[1] = key[1]
        q[2] = key[2]

    if key[0] == '':
        q[0] = '*'
    else:
        q[0] = key[0]

    query = f".//{q[0]}{q[1]}{q[2]}"
    return query


def remove_element(root, query):
    element_ids = root.findall(query)

    if isinstance(element_ids, list):
        for element in element_ids:
            for parent in root.iter():
                for child in list(parent):
                    if child is element:
                        parent.remove(child)

    else:
        for element in element_ids:
            root.remove(element)

    return root


def add_empiricaltree(root, tree_path):
    new_sub1 = ET.Element('empiricalTreeDistributionModel')
    new_sub1.set("id", "treeModel")
    new_sub1.set("fileName", tree_path)
    new_sub1a = ET.SubElement(new_sub1, 'taxa')
    new_sub1a.set("idref", "taxa")

    # after taxa but before 'strictClockBranchRates' = additional flexibility required
    position = root.find('taxa')
    root.insert(list(root).index(position) + 1, new_sub1)

    new_sub2 = ET.Element('statistic')
    new_sub2.set("id", "treeModel.currentTree")
    new_sub2.set("name", 'Current Tree')
    new_sub2a = ET.SubElement(new_sub2, 'empiricalTreeDistributionModel')
    new_sub2a.set("idref", "treeModel")
    position = root.find('taxa')
    root.insert(list(root).index(position) + 2, new_sub2)

    target_element = root.find(".//operators")
    new_sub3 = ET.SubElement(target_element, 'empiricalTreeDistributionOperator')
    new_sub3.set("weight", "3")
    new_sub3a = ET.SubElement(new_sub3, 'empiricalTreeDistributionModel')
    new_sub3a.set("idref", "treeModel")

    return root


def main():
    inputs = parse_args()
    xml = ET.parse(inputs.xml_path)
    root = xml.getroot()
    file_path = re.sub('.xml', '_empiricaltree.xml', inputs.xml_path)

    for key in keys:
        query = parse_query(key)
        root = remove_element(root, query)

    root = add_empiricaltree(root, inputs.tree_path)

    xml.write(file_path)



if __name__ == "__main__":
    main()
