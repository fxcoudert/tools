#!/usr/bin/env python3

import json
import os
import sys

import requests


class ChemrxivClient:
    """Handle figshare API requests, using access token."""

    base = 'https://api.figshare.com/v2'
    pagesize = 100

    def __init__(self, token):
        """Initialiase the object and check access to the API."""

        self.token = token
        self.headers = {'Authorization': 'token ' + self.token}

        r = requests.get(f'{self.base}/account', headers=self.headers)
        r.raise_for_status()

    def request(self, url, method, params):
        """Send a figshare API request."""
        if method.casefold() == 'get':
            return requests.get(url, headers=self.headers, params=params)
        elif method.casefold() == 'post':
            return requests.post(url, headers=self.headers, json=params)
        else:
            raise Exception(f'Unknown method for query: {method}')

    def query(self, query, method='get', params=None):
        """Perform a direct query."""
        r = self.request(f'{self.base}/{query}', method, params)
        r.raise_for_status()
        return r.json()

    def query_generator(self, query, method='get', params=None):
        """Query for a list of items, with paging and return a generator."""
        if params is None:
            params = {}
        n = 0
        while True:
            params.update({'limit': self.pagesize, 'offset': n})
            r = self.request(f'{self.base}/{query}', method, params)
            r.raise_for_status()
            r = r.json()

            # Special case if a single item, not a list, was returned
            if not isinstance(r, list):
                yield r
                return

            # If we have no more results, bail out
            if len(r) == 0:
                return

            yield from r
            n += self.pagesize

    def query_list(self, *args, **kwargs):
        """Query of a list of item, handling paging internally, returning a list. 

        May take a long time to return."""
        return list(self.query_generator(*args, **kwargs))

    def all_preprints(self):
        """Return a generator to all the chemRxiv preprints."""
        return self.query_generator('articles?institution=259')

    def preprint(self, identifier):
        """Information on a given preprint."""
        return self.query(f'articles/{identifier}')

    def author(self, identifier):
        """Information on a given preprint."""
        return self.query(f'account/authors/{identifier}')

    @staticmethod
    def custom_fields_as_dict(doc):
        """Retrieve chemRxiv custom fields as a dictionary."""
        return {i['name']: i['value'] for i in doc['custom_fields']}

    def search_authors(self, criteria):
        """Search for authors."""
        return self.query('account/authors/search', method='POST', params=criteria)

    def search_preprints(self, criteria):
        """Search for preprints."""
        p = {**criteria, 'institution': 259}
        return self.query_list('articles/search', method='POST', params=p)


# Below, we demonstrate how this works
def main():
    # We need a Figshare API token, see https://docs.figshare.com
    # You can insert it below, or store it as text in ~/.figshare_token
    token = 'invalid'
    try:
        with open(os.path.expanduser('~/.figshare_token'), 'r') as file:
            token = file.read().strip()
    except IOError:
        pass

    # Connect to Figshare
    try:
        api = ChemrxivClient(token)
    except requests.exceptions.HTTPError as e:
        print(f'Authentication did not succeed. Token was: {token}')
        print(f'Error: {e}')
        sys.exit(1)

    # Find details of a given preprint
    doc = api.preprint(5349151)
    print(f'Title for preprint {doc["id"]} is: {doc["title"]}')
    print('Full data available:')
    print(json.dumps(doc, indent=4, sort_keys=True))
    print('\n')

    # Find author details
    author_id = doc['authors'][0]['id']
    author = api.author(author_id)
    print(f'First author is {author["full_name"]}')
    print('Full data available:')
    print(json.dumps(author, indent=4, sort_keys=True))
    print('\n')

    # We can iterate over all preprints
    for doc in api.all_preprints():
        if 'zeolite' in doc['title'].lower():
            # Get more details
            print(f'Found a zeolite paper! Title is: {doc["title"]}')
            details = api.preprint(doc['id'])
            authors = {x['full_name'] for x in details['authors']}
            print(f'Its authors are: {", ".join(authors)}')

            # Iterating over preprints, we stop now
            break
    print('\n')

    # We can also search for authors
    print('Searching for author with given ORCID: 0000-0002-8277-4434')
    res = api.search_authors({'orcid': '0000-0002-8277-4434'})
    print(f'Found {len(res)}:')
    for i in res:
        print(f'Name: {i["full_name"]}, job: {i["job_title"]}')
    print('\n')

    print('Searching for authors named Jennifer:')
    res = api.search_authors({'search_for': 'Jennifer'})
    print(f'Found {len(res)}:')
    for i in res:
        print(f'  - {i["full_name"]}')
    print('\n')

    # And of course, search for preprints
    res = api.search_preprints({'search_for': 'zeolite'})
    print(f'Found {len(res)} papers about “zeolite”, here are the first five:')
    for i in res[:5]:
        print(f'  - {i["title"]}')

    sys.exit(0)


if __name__ == '__main__':
    main()
