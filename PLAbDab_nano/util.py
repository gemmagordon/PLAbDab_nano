import pandas as pd


def add_url_to_data(df):
    links = []

    for entry in df.iloc:

        if "patent US" in entry.definition:
            patent_id = entry.definition.split()[-1]
            link = "https://patents.google.com/patent/US{}/en".format(patent_id)
        # NOTE changed 'pairing' col to 'source' 
        elif entry.source == "SAbDab":
            link = "https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/structureviewer/?pdb={}".format(entry.ID[:4].lower())
        elif entry.source == "TheraSAbDab":
            link = "https://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/search/?therapeutic={}".format(entry.ID.replace("_2",""))
        
        elif entry.source == "GenBank":
            link = "https://www.ncbi.nlm.nih.gov/protein/{}".format(entry.ID)

        links.append(link)

    df = df.copy()
    df["url"] = links
    return df