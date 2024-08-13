import socket
import time
s = socket.socket(socket.AF_INET,socket.SOCK_DGRAM,0)
sip = "192.168.3.121"
sport = 27181
s.connect((sip,sport))
dest = socket.gethostname()
print dest
#c, addr = s.accept()

#print addr
#s.bind(("192.168.3.121",27181))
s.settimeout(.2)
'''
while(2 >0):
	sdata = raw_input ('>')
	sdata = sdata.decode('hex')
	s.sendto(sdata,(sip,sport))
	try:
		data,addr = s.recvfrom(1280)
		print ('>'),data.encode('hex')
	except socket.timeout:
		print ('No answer')
'''
i=0
command_write_1 = hex(0x81c20820)
data_write_lb = hex(0x00)
data_write_hb = hex(0x00)
max_low = hex(0xff)
command_read = hex(0x01420820)
#command_read = command_read.decode('hex')
while(int(data_write_hb,16)!=int(max_low,16)):
	#a=str(command_read)
	#a = command_read.decode('hex')
	#print hex(command_read)

	s.sendto(command_read[2:].zfill(8).decode('hex'),(sip,sport))
	try:
		data,addr = s.recvfrom(1280)
		print ('Received: '),data.encode('hex')
	except socket.timeout:
		print ('No answer')

	if(int(data_write_lb,16)+1>int(max_low,16)):
		data_write_hb = hex(int(data_write_hb,16)+1)
		data_write_lb = hex(0x00)
	
	packet = 0
	data_write_lb = hex(int(data_write_lb,16) + 1)
	packet = hex(int(command_write_1,16)<<8 | int(data_write_lb,16)) 
	packet = hex(int(packet,16)<<8 | int(data_write_hb,16))
	packet[2:].decode('hex')
	time.sleep(0.1)

	s.sendto(packet[2:].zfill(12).decode('hex'),(sip,sport))
	try:
		data,addr = s.recvfrom(1280)
	except socket.timeout:
		print ('Wrote: '),packet[-4:]
		

	time.sleep(0.1)
	i=i+1