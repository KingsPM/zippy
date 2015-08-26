# primula
Primer database and design

## Description
This program integrates a simple SQLite primer database and design tool.

## Usage
To design primers and query existing
> `primula.py <VCF> --design`

add primers to database
> `primula.py <PRIMERS> --add`

retrieve (and design) primer for a single location
> `primula.py -q chr1:202030 --design`

