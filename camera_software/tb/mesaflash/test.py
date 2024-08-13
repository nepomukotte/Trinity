import socket
s = socket.socket(socket.AF_INET,socket.SOCK_DGRAM,0)
sip = "192.168.1.121"
sport = 27181
s.settimeout(.2)
while(2 >0):
	sdata = raw_input ('>')
	sdata = sdata.decode('hex')
	s.sendto(sdata,(sip,sport))
	try:
		data,addr = s.recvfrom(1280)
		print ('>'),data.encode('hex')
	except socket.timeout:
		print ('No answer')
