cd /local/httpd/custom

echo --- before stop ---
ps -ef | grep httpd
sudo ./stop_apache.sh
sudo ./stop_apache.sh
echo --- after stop ---
ps -ef | grep httpd
sleep 3;
sudo ./start_apache.sh
echo --- after start ---
ps -ef | grep httpd