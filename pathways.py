#!/usr/bin/python

from intermine.webservice import Service
from intermine.webservice import ServiceError
from collections import defaultdict
from itertools import groupby, ifilter
import org_util

class PathwayDemo(object):
	service_urls = {'Drosophila melanogaster': "http://www.flymine.org/query/service",
			'Saccharomyces cerevisiae': "http://yeastmine.yeastgenome.org/yeastmine/service",
			'Rattus norvegicus': "http://ratmine.mcw.edu/ratmine/service",
			'Mus musculus': "http://metabolicmine.org/test/service",
			'Homo sapiens': "http://metabolicmine.org/test/service",
	}

	def __init__(self):
		self.services = {}
		for (name, service_url) in self.service_urls.items():
			try:
				self.services[name] = Service(service_url)
				print "Connected to %s" % service_url
			except ServiceError:
				print "Failed to initialise: %s" % service_url
		 

	# returns a list of lists
	def find_gene(self, symbol, org_name):
		service = self.services[org_name]

		query = service.new_query("Gene")\	
			       .select("symbol", "primaryIdentifier", "name", "organism.name")\
		               .where("Gene", "LOOKUP", symbol)\
		               .where("organism.name", "=", org_name)

		return [row.to_l() for row in query.rows()]
				

	def get_homologs_for_gene(self, symbol, org_name):
		# always use FlyMine for querying homologs
		
		h_sym = "homologues.homologue.symbol"
		h_org = "homologues.homologue.organism.name"
		h_ds  = "homologues.dataSets.name"
		
		service = self.services["Drosophila melanogaster"]
    	
		query = service.new_query("Gene")\
		               .select(h_org, h_sym, h_ds)\
		               .where(h_org, "ONE OF", org_util.get_names())\
		               .where("organism.name", "=", org_name)\
		               .where("symbol", "=", symbol)\
		               .where(h_sym, "IS NOT NULL")

		homologs = defaultdict(dict)		
		
		for org, g1 in groupby(query.rows(), lambda x: x[h_org]):
			for sym, g2 in groupby(g1, lambda x: x[h_sym]):
				homlogues[org][sym] = [r[h_ds] for r in g2]
		return homologs
 
	def get_pathways(self, symbol, org_name):
		service = self.services[org_name]
		sym = "symbol"
		pw = "pathways.name"
		ds = 'Gene.pathways.dataSets.name'
		org = "organism.name"
		query = service.new_query("Gene").select(sym, org, pw)
	
		# YeastMine doesn't have pathway.dataSets, check model first
		if self.is_path_in_model(org_name, ds):
			query.add_view(ds)
			query.add_join('Gene.pathways.dataSets', 'OUTER')

		query.add_sort_order(pw, 'asc')
		query.add_constraint(sym, '=', symbol)
		query.add_constraint(org, "=", org_name)
	
	        # Return a list of triples
		return [(r[pw] + " (%s)" % r[ds] if len(r) == 3 else r[pw], r[org], r[sym]) for r in query.rows()]

	def is_path_in_model(self, org, path):
		service = self.services[org]	
		try:
			service.model.validate_path(path)
		except:
			return False
		return True

	def strip_suffix(self, dataset):
		if dataset.endswith('data set'):
			return dataset[0:dataset.find('data set')]
