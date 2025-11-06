import sys
from pathlib import Path

import requests


class NcbiApi:
    def __init__(self, key):
        self.ncbi_api_uri = "https://api.ncbi.nlm.nih.gov/datasets/v2"
        self.ncbi_api_key = key

    def get_request(self, command: str) -> requests.Response:
        """
        Returns a requests response object for a given NCBI api query (command).

        args:
            command -> str: request to append to NCBI API URI
        """
        if self.ncbi_api_key != "":
            headers = {"Accept": "application/json", "api-key": f"{self.ncbi_api_key}"}
        else:
            headers = {"Accept": "application/json"}

        response = requests.get(self.ncbi_api_uri + command, headers=headers)

        if response.status_code != 200:
            raise Exception(
                f"Cannot connect to NCBI (status code '{str(response.status_code)}')'"
            )

        return response

    def get_assemblies_for_taxon(
        self, taxon: str, filters_reference_only: bool = False, page_size: int = 1000
    ):
        """
        For a given taxon name or taxid, a get a list describing all available
        assemblies for a taxon.

        args:
            taxon -> string: Taxon name or taxid to query
            filters_reference_only -> bool: Return reference assemblies only
            page_size -> int: Number of assembies to return per page (max 1000)
        """
        query = f"/genome/taxon/{taxon}/dataset_report"
        arg_filter_refs_only = f"?filters_reference_only={filters_reference_only}"
        arg_page_size = f"&page_size={page_size}"
        r = self.get_request(query + arg_filter_refs_only + arg_page_size)
        r_res = r.json()

        assemblies = []

        if r_res.get("reports"):
            assemblies = [] + r_res.get("reports")

        while r_res.get("next_page_token", "") != "":
            page_token_arg = f"&page_token={r_res.get('next_page_token')}"
            r = self.get_request(
                query + arg_filter_refs_only + arg_page_size + page_token_arg
            )
            if r_res.get("reports", None):
                assemblies = assemblies + r_res.get("reports")

        return assemblies

    def assembly_count_for_taxon(
        self, taxon: str, filters_assembly_source: str = "all"
    ) -> int:
        """
        For a given taxon name or taxid, return the number of assemblies available
        on NCBI.

        args:
            taxon -> string: Taxon name or taxid to query
            filters_assembly_source -> string: "refseq", "genbank", or "all"
        """
        query = f"/genome/taxon/{taxon}/dataset_report"
        filter_source = f"?filters_assembly_source={filters_assembly_source}"
        page_size = "&page_size=1"
        taxon = self.get_request(query + filter_source + page_size)

        return taxon.json().get("total_count", 0)

    def download_genomes(self, accessions: list, outfile: Path) -> Path:
        if self.ncbi_api_key != "":
            headers = {
                "Accept": "application/zip",
                "api-key": f"{self.ncbi_api_key}",
                "content-type": "application/json",
            }
        else:
            headers = {"Accept": "application/zip", "content-type": "application/json"}

        try:
            r = requests.post(
                self.ncbi_api_uri + "/genome/download",
                headers=headers,
                json={
                    "accessions": accessions,
                    "include_annotation_type": ["GENOME_FASTA"],
                    "hydrated": "FULLY_HYDRATED",
                },
                stream=True,
            )
            r.raise_for_status()

            downloaded = 0
            with open(outfile, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        print(
                            f"\rDownloading: {downloaded / 1024 / 1024:.1f}MB", end=""
                        )

            return outfile

        except requests.HTTPError as e:
            print(f"Failed to download genomes for accessions: {accessions}")
            sys.exit(f"Reason: {e}")
