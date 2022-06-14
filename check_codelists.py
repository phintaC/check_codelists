#################################################################
#   Author:     Josephine Strange (JS)                          #
#   Program:    check_codelists.py                              #
#   Date:       2022-06-14                                      #
#   Version:    1.0                                             #
#   Notes:      Does not recognise integer values in pandas     #
#               DataFrames, converting them to floats, why      #
#               integer-based codelists will hit                #
#################################################################
import re, sys, os
import xml.etree.ElementTree as ET 
import pyreadstat as prs # Faster than pandas in reading SAS files

ipath = input("Please enter root of study: ")
isdtm = ipath + r"\\Data\\SDTM"
ixml = isdtm + r"\\xpt\\define.xml"

# Parse the XML file
print("Parsing define.xml...")
tree = ET.parse(ixml)

################################## Fetch ItemDef (CodeList reference) ###################################
items,codelists = dict(),dict()
for node in tree.findall('.//{http://www.cdisc.org/ns/odm/v1.3}ItemDef'):
    itemOID = node.attrib['OID']

    for child in node.iter():
        if re.search(r'CodeListRef',child.tag):
            if child.attrib['CodeListOID'] not in {"CL.DRUGDICT", "CL.MEDDRA"}:
                codelistOID = child.attrib['CodeListOID']
            else: 
                continue

            if itemOID not in items:
                items[itemOID] = [codelistOID,None,None]
            else:
                items[itemOID] += [codelistOID,None,None]

            if codelistOID not in codelists: 
                codelists[codelistOID] = [itemOID]
            else:
                codelists[codelistOID] += [itemOID]

############################# Fetch WhereClauses ####################################
whereclauses = dict()
for node in tree.findall('.//{http://www.cdisc.org/ns/def/v2.0}ValueListDef'):
    # get WhereClause names
    for child in node.iter():
        # Get the item reference 
 
        if re.search(r"ItemRef",child.tag):
            itemOID = child.attrib["ItemOID"]

        elif re.search(r'WhereClauseRef',child.tag):
            whereClauseOID = child.attrib['WhereClauseOID']
            if whereClauseOID not in whereclauses: 
                whereclauses[whereClauseOID] = [itemOID]
            else:
                whereclauses[whereClauseOID] += [itemOID]

################################ WhereClause conditions #################################
whereclauses_na = []
for node in tree.findall('.//{http://www.cdisc.org/ns/def/v2.0}WhereClauseDef'):
    n = 0
    whereClauseOID = node.attrib['OID']

    # Fetch the conditions - used for eval() statement
    for child in node.iter():
        if re.search(r"RangeCheck",child.tag):

            conditional_variable = ".".join(child.attrib['{http://www.cdisc.org/ns/def/v2.0}ItemOID'].split(".")[2:])
            comparator = child.attrib["Comparator"].replace("EQ","==").replace("NE","!=").replace("NOTIN", "not in").lower()
            checkValue = [grandchild.text for grandchild in child.iter() if re.search(r"CheckValue",grandchild.tag)]

            if "," in checkValue[0]: 
                checkValue = checkValue[0].replace(" ","").split(",")

            # Make conditions python readable; distinguish between strings and lists 
            if len(checkValue) == 1:
                value = "'" + "".join(checkValue) + "'"
            else:
                value = ".isin(['" +  "','".join(checkValue) + "'])"

            # Make python compilable and concatenate conditions if there is an "and" statement
            if n == 0:
                if comparator == "not in":
                    condition = "(~df['" + conditional_variable +"']" + value + ")"
                elif comparator == "in":
                    condition = "(df['" + conditional_variable +"']" + value + ")"
                else: 
                    condition = "(df['" + conditional_variable + "'] " + comparator + " " + value + ")"

            else:
                if comparator == "not in":
                    condition += " & (~df['" + conditional_variable + "']" + value + ")"
                elif comparator == "in":
                    condition += " & (df['" + conditional_variable + "']" + value + ")"
                else: 
                    condition += " & (df['" + conditional_variable + "'] " + comparator + " " + value + ")"
                
            n += 1

    # # Add conditions to items dictionary
    if whereClauseOID in whereclauses: 
        itemOIDs = whereclauses[whereClauseOID]

        for itemOID in itemOIDs:
            if itemOID in items:
                items[itemOID][1] = condition.replace(" in ","")
        
        whereclauses.pop(whereClauseOID)
    else:
        whereclauses_na += [whereClauseOID] 
        print("WARNING: WhereClause not accounted for: " + whereClauseOID)
        continue

################################ Fetch CodeLists #################################
for node in tree.findall('.//{http://www.cdisc.org/ns/odm/v1.3}CodeList'):
    if node.attrib['OID'] not in {"CL.DRUGDICT", "CL.MEDDRA"}:
        codelistOID = node.attrib['OID']
    else:
        continue

    codes = []
    for child in node.iter():
        if re.search(r'CodeListItem',child.tag):
            codes += [child.attrib['CodedValue']]
        
    for itemOID in codelists[codelistOID]:
            items[itemOID][2] = codes
            items[itemOID].pop(0) # Removes the codelist reference

################################ Fetch domain and variables #################################
for key in list(items.keys()): 
    domain_variable = ".".join(key.split(".")[1:3])

    # Create a list of lists per domain.variable; 
    if domain_variable not in items:
        items[domain_variable] = [items[key]]
    else:
        items[domain_variable] += [items[key]]

    items.pop(key)

# Discard the obsolete dictionaries 
whereclauses.clear()
codelists.clear()

# Create set of domain : variables to be checked
domainVariables = dict()
for key in list(items.keys()):
    if key.split(".")[0] not in domainVariables:
        domainVariables[key.split(".")[0]] = {key.split(".")[1]}
    else:
        domainVariables[key.split(".")[0]].add(key.split(".")[1])

################################ Start looking in datasets #################################
print("Checking codelists...")

try: 
    files = [f for f in os.listdir(isdtm) if os.path.isfile(isdtm + r"\\" + f) and f.split(".")[0].upper() in {key.split(".")[0] for key in items.keys()}]

except WindowsError as err:
    print("Error: " + err)
    sys.exit(1)

missing_items,missing_values = dict(), dict()

for f in files:
    df,meta = prs.read_sas7bdat(isdtm + r"\\" + f)

    variables = domainVariables[f.split(".")[0].upper()]

    # Check if variables listed in define-XML exists in datasets, or unique set of values from dataset exists in define-XML
    for v in variables:
        domain_variable = f.split(".")[0].upper() + "." + v
        
        if domain_variable in items: 
            for item in items[domain_variable]:
                define_codelist = item[1]
                data_codelist = set(df[v])

                # Unconditional: 
                if item[0] == None:
                    missing_item = ["", list( set(item[1]) - set(df[v].unique()) )] # Missing in dataset
                    missing_value = ["", list( set(df[v].unique()) - set(item[1]) )] # Missing in define-XML

                # Conditional
                else: 
                    subset = df.loc[(eval(item[0]))] 
                    missing_item = [item[0], list( set(item[1]) - set(subset[v].unique()) )] # Missing in dataset
                    missing_value = [item[0], list( set(subset[v].unique()) - set(item[1]) )] # Missing in define-XML

                # if "" in missing_item[1]:
                #     missing_item[1].remove("")
                if "" in missing_value[1]:
                    missing_value[1].remove("")

                # Add the missing items/values to the dictionaries of missing items/values only if list of missing items/values is not empty
                if len(missing_value[1]) > 0:
                    if domain_variable not in missing_values:
                        missing_values[domain_variable] = [missing_value]
                    else:
                        missing_values[domain_variable] += [missing_value]

                if len(missing_item[1]) > 0:
                    if domain_variable not in missing_items:
                        missing_items[domain_variable] = [missing_item]
                    else:
                        missing_items[domain_variable] += [missing_item]


with open(ipath + r"\\codelist_check.txt","w+") as outfile:
    outfile.write("Please note: integer-based variables will be read as float values and therefore, these will hit as missing in dataset/define-XML")
    outfile.write("\n\t\t\"nan\" is printed to CODELIST VALUES when codelist from define-XML contains character values, while the variable is numeric")
    outfile.write("\n\t\tThe way WHERECLAUSE CONDITIONs are represented is to make the condition python readable")
    
    if len(whereclauses_na) > 0:
        outfile.write("\n\n" + "#"*25 + " WhereClauses not accounted for " + "#"*25 + "\n")
        outfile.write("\n".join(whereclauses_na))
        outfile.write("\n##################################################################################")

    if len(missing_items) > 0 or len(missing_values) > 0:
        outfile.write("\n\nDATASET.VARIABLE\tWHERECLAUSE CONDITION\nCODELIST CODES")

    if len(missing_items) > 0:
        outfile.write("\n" + "#"*25 + " Missing in dataset " + "#"*25)
        for key in missing_items:
            for item in missing_items[key]:
                outfile.write("\n" + key + ":\t" + item[0] + "\n" + str(item[1])[1:-1] + "\n")

    if len(missing_values) > 0:
        outfile.write("\n\n" + "#"*25 + " Missing in define-XML " + "#"*25)
        for key in missing_values:
            for item in missing_values[key]:
                outfile.write("\n" + key + ":\t" + item[0] + "\n" + str(item[1])[1:-1] + "\n")