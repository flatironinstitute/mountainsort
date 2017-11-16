Plans for cordion authentication server
=======================================

As discussed with Witold on 10/11/2017, here's the current plan for Cordion authentication server for remote processing using MountainLab:

* I am the admin of river (processing server). 
* I register with the cordion web site, asking for approval. my application is accepted. hurray!
* I am given a private key associated with river. Only cordion and I know it.
* Cordion adds 'river'/'river-public-key' to its public table of registered servers 
* Now I configure my server pasting 'river'/'river-private-key' into a configuration file and launch kulelepoller which tries to connect to kulele
* Kulele gets the connection request which is signed using 'river-private-key', and then it verifies that it was properly signed by comparing with cordion's public table. Kulele accepts the connection.
* The browser (client) requests from cordion authorization for some@user on river (proving that some@user has authenticated with a service like google).
* Cordion returns the authorization token signed using 'river-private-key'
* That token is sent along with all job requests to kulele
* Kulele verifies the signature is valid for river (again by consulting cordion's public table) 
* Kulele checks the request is consistent with the authorization and forwards the request to river
* River trusts kulele and executes the processing job