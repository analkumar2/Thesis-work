#exec(open('celllist.py').read())

#Execute for all the html files

import html2text
import urllib3
import os
import re
from bs4 import BeautifulSoup
import requests

# For stored html file
file = open('https _channelpedia.epfl.ch_expdata3.html')
astext = file.read()

# #For online html
# http = urllib3.PoolManager()
# url = 'https://channelpedia.epfl.ch/expdata'
# response = http.request('GET', url)
# soup = BeautifulSoup(response.data)
# astext = str(soup)

listcellid = re.findall(r'Cell ID \d+',astext)
listid = []
for cellid in listcellid:
    listid.append(re.findall(r'\d+',cellid)[0])
print(listid)

file.close()

urllink_list = [f'https://channelpedia.epfl.ch/expdata/download/{id}/rCell/nwb' for id in listid]
print(urllink_list)

urllink_str = ''
for url in urllink_list:
    urllink_str = urllink_str+' \n'+url

file = open('urllink_str.txt', 'a')
file.write(urllink_str)
file.close()

#Now use the urllink_str.txt file with any in-browser batch downloader
