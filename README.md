## prostPostprocess

A R package providing functions to enable differential analysis of Prost! output

### Installation

The package depends of DESeq2. DESeq2 can be installed from Bioconductor. Other depedencies can be automatically installed by the package manager from the CRAN repository.

The package can be installed from the github repository if you have access to it : 

```
devtools::install_github(repo = "https://github.com/INRAE-LPGP/prostPostprocess", auth_token = "your_PAT")
```

To generate a personnal access token (PAT), go to https://github.com/settings/tokens and supply to `auth_token`
