# dependencies
# pip install pandas aiohttp tqdm
import asyncio
import urllib.parse
import aiohttp
import pandas as pd
from tqdm.asyncio import tqdm


async def fetch_name_from_smiles(session, smiles, semaphore):
    """Worker function that queries PubChem for a single SMILES string

    while respecting a concurrency limit (semaphore).
    """
    if pd.isna(smiles) or str(smiles).strip() == "":
        return "Missing SMILES", "None"

    # URL encode special characters like # and @ in SMILES strings
    safe_smiles = urllib.parse.quote(str(smiles).strip())
    prop_url = f"https://nih.gov{safe_smiles}/property/IUPACName/json"

    # Control concurrency using the semaphore
    async with semaphore:
        try:
            # Step 1: Look up CID and IUPAC Name
            async with session.get(prop_url, timeout=10) as response:
                if response.status == 200:
                    data = await response.json()
                    prop_data = data["PropertyTable"]["Properties"][0]
                    cid = prop_data.get("CID")
                    iupac_name = prop_data.get("IUPACName")

                    if cid:
                        # Step 2: Look up Synonyms using the discovered CID
                        syn_url = f"https://nih.gov{cid}/synonyms/json"
                        async with session.get(syn_url, timeout=10) as syn_resp:
                            if syn_resp.status == 200:
                                syn_data = await syn_resp.json()
                                synonyms = syn_data["InformationList"][
                                    "Information"
                                ][0].get("Synonym", [])
                                if synonyms:
                                    # Fallback achieved: Found common name!
                                    return synonyms[0], "common"

                    # Fallback achieved: Common name missing, using IUPAC
                    if iupac_name:
                        return iupac_name, "iupac"
                    return "Name Unavailable", "none"

                elif response.status == 404:
                    return "Not Found in Database", "none"
                else:
                    return f"API Error ({response.status})", "none"

        except Exception as e:
            return f"Error: {str(e)}", "none"


async def main_async(input_csv, output_csv, smiles_column):
    # Load your CSV data
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)

    if smiles_column not in df.columns:
        raise ValueError(f"Column '{smiles_column}' not found in the CSV!")

    # PubChem rules: Max 5 requests/sec per IP, but their servers handle bursts well
    # A semaphore of 5 limits us to 5 concurrent active connections at any exact millisecond
    sem = asyncio.Semaphore(5)

    # Establish the persistent HTTP session
    async with aiohttp.ClientSession() as session:
        tasks = [
            fetch_name_from_smiles(session, smiles, sem)
            for smiles in df[smiles_column]
        ]

        print(f"Starting async processing for {len(tasks)} molecules...")

        # Gather tasks together and show an interactive progress bar
        results = await tqdm.gather(*tasks, desc="Resolving Names")

    # Unpack the tuple results into your dataframe columns
    df["resolved_name"] = [r[0] for r in results]
    df["name_type"] = [r[1] for r in results]

    # Save data
    df.to_csv(output_csv, index=False)
    print(f"\nSuccess! Output saved to: {output_csv}")


if __name__ == "__main__":
    # Local Configuration
    INPUT_FILE = "my_molecules.csv"  # Path to your input file
    OUTPUT_FILE = "resolved_molecules.csv"  # Path to save results
    SMILES_COL_NAME = "smiles"  # Change this to match your CSV column header

    # Run the asynchronous loop
    asyncio.run(
        main_async(INPUT_FILE, OUTPUT_FILE, smiles_column=SMILES_COL_NAME)
    )

