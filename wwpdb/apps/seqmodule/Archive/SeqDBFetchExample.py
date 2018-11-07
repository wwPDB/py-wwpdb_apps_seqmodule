use FetchUnpXml.py:

obj = FetchUnpXml(id)
dict = obj.ParseUnpXmlData()

where id is 'PDGFB_HUMAN' or 'P01127'. using name or accession id gets same result.


dict has key value pair for the following keys:

'db_code'
'db_accession'
'sequence'
'ec'
'keyword'
'name'
'synonyms'
'gene'
'source_scientific'
'source_common'
'taxonomy_id'
'comments'

similar for genebank RNA sequence, use FetchNcbiXml.py

obj = FetchNcbiXml(id)
dict = obj.ParseNcbiXmlData()

where id should be gi number.

dict has key value pair for the following keys:

'sequence'
'source_scientific'
'taxonomy_id'
