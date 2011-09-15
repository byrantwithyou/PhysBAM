# Sim Configuration

import os.path

# session directory (location independent) never moved
session_directory="/n/curvature/data/psim/data/sessions"
# SSL 
keys_directory="/n/curvature/data/psim/keys"
ca_certificate_file=os.path.join(keys_directory,"CA.cert")

server_private_key_file=os.path.join(keys_directory,"server.pkey")
server_certificate_file=os.path.join(keys_directory,"server.cert")

client_private_key_file=os.path.join(keys_directory,"client.pkey")
client_certificate_file=os.path.join(keys_directory,"client.cert")
# hosts file in form specified in SERVER
hosts_server="curvature.stanford.edu"
hosts_port=7777
# port to listen upon
server_port=8888
