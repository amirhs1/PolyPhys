#!/bin/bash
# 1.Do this on your machine
ssh-keygen # here you need to confimr key generation and setup a password for your rsa key
cd ssh; ls
# 2.know you login to your compute canada account
sftp username@graham.computecanada.ca
# 3a.athe following are done if you have not have already an rsa key.
mkdir .ssh
chmod 700 .ssh
cd .ssh
# 3b.do this for your new key
cd .ssh
# 4. final step
put id_rsa.pub authorized_keys
chmod 644 authorized_keys
exit
# 5. now test it
ssh username@graham.computecanada.ca

