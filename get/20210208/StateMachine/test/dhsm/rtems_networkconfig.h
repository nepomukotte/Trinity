/*
 * P. Kestener (8/06/2008) 
 *
 * principales modif pou le BSP Virtex:
 * - modifier ethernet_address
 * - rtems_bsdnet_config
 * - netdriver_config en utilisant le
 *   le message suivant
 *   http://osdir.com/ml/os.rtems.user/2007-05/msg00124.html
 *
 */

/*
 * Network configuration -- VIRTEX BSP
 *
 * See one of the other networkconfig.h files for an
 * example of a system that includes a real NIC and
 * the loopback interface.
 * 
 ************************************************************
 * EDIT THIS FILE TO REFLECT YOUR NETWORK CONFIGURATION     *
 * BEFORE RUNNING ANY RTEMS PROGRAMS WHICH USE THE NETWORK! * 
 ************************************************************
 *
 *  $Id: networkconfig.h,v 1.13 2008/02/05 21:49:36 joel Exp $
 */

#ifndef _RTEMS_NETWORKCONFIG_H_
#define _RTEMS_NETWORKCONFIG_H_

/*
 * default behaviour is to get boot parameter (IP adress) from a
 * bootp/dhcp server.
 * Comment the following line for a manual network configuration and
 * edit structure netdriver_config, defined right below.
 * Take care to provide the MAC address of the RTEMS target.
 */
#define RTEMS_USE_BOOTP

/*
 *  The following will normally be set by the BSP if it supports
 *  a single network device driver.  In the event, it supports
 *  multiple network device drivers, then the user's default
 *  network device driver will have to be selected by a BSP
 *  specific mechanism.
 */

#if defined(leon2)
#define RTEMS_BSP_NETWORK_DRIVER_NAME   "smc_91111"
#define RTEMS_BSP_NETWORK_DRIVER_ATTACH RTEMS_BSP_NETWORK_DRIVER_ATTACH_SMC91111
#endif // leon2

#if defined(leon3)
#define RTEMS_BSP_NETWORK_DRIVER_NAME   RTEMS_BSP_NETWORK_DRIVER_NAME_GRETH
#define RTEMS_BSP_NETWORK_DRIVER_ATTACH RTEMS_BSP_NETWORK_DRIVER_ATTACH_GRETH
#warning "leon3"
#endif // leon3

#if defined(pc386)
/*
 * Network configuration -- QEMU using DHCP
 */
#define RTEMS_BSP_NETWORK_DRIVER_NAME    "ne1"
#define RTEMS_BSP_NETWORK_DRIVER_ATTACH  rtems_ne_driver_attach 
#undef RTEMS_USE_BOOTP
#endif // pc386

/*
 * check network driver is defined
 */
#ifndef RTEMS_BSP_NETWORK_DRIVER_NAME
#warning "RTEMS_BSP_NETWORK_DRIVER_NAME is not defined"
#define RTEMS_BSP_NETWORK_DRIVER_NAME "no_network1"
#endif

#ifndef RTEMS_BSP_NETWORK_DRIVER_ATTACH
#warning "RTEMS_BSP_NETWORK_DRIVER_ATTACH is not defined"
#define RTEMS_BSP_NETWORK_DRIVER_ATTACH 0
#endif

#include <bsp.h>

/*
 * Define RTEMS_SET_ETHERNET_ADDRESS if you want to specify the
 * Ethernet address here.  If RTEMS_SET_ETHERNET_ADDRESS is not
 * defined the driver will choose an address.
 */
#define RTEMS_SET_ETHERNET_ADDRESS
#if (defined (RTEMS_SET_ETHERNET_ADDRESS))
/*
 * the following mac address is used for the leon2/at697 Gaisler's
 * board (SVOM/UTS testbench).
 */
static char ethernet_address[6] = { 0x12, 0xC2, 0x6F, 0x95, 0x04, 0xD9 };
#endif


#ifdef RTEMS_USE_LOOPBACK 
/*
 * Loopback interface
 */
extern int rtems_bsdnet_loopattach();
static struct rtems_bsdnet_ifconfig loopback_config = {
  "lo0",                    /* name */
  rtems_bsdnet_loopattach,  /* attach function */
  NULL,                     /* No more interfaces */
  "127.0.0.1",              /* IP address */
  "255.0.0.0",              /* IP net mask */
  NULL,                     /* Driver supplies hardware address */
  0,                        /* Use default driver parameters */
  0,                        /* default efficiency multiplier */
  0,                        /* default udp TX socket buffer size */
  0,                        /* default udp RX socket buffer size */
  0,                        /* default tcp TX socket buffer size */
  0,                        /* default tcp RX socket buffer size */
};
#endif /* RTEMS_USE_LOOPBACK */

/*
 * Default network interface
 */
static struct rtems_bsdnet_ifconfig netdriver_config = {
  RTEMS_BSP_NETWORK_DRIVER_NAME,		/* name */
  RTEMS_BSP_NETWORK_DRIVER_ATTACH,	/* attach function */

#ifdef RTEMS_USE_LOOPBACK 
  &loopback_config,		/* link to next interface */
#else
  NULL,				/* No more interfaces */
#endif /* RTEMS_USE_LOOPBACK */

#if (defined (RTEMS_USE_BOOTP))
  NULL,				/* BOOTP supplies IP address */
  NULL,				/* BOOTP supplies IP net mask */
#else
#  if defined(pc386)
  "10.0.2.5",                   /* IP address of target inside qemu */
#  else
  "132.166.9.214",		/* IP address */
#  endif
  "255.255.0.0",		/* IP net mask */
#endif /* !RTEMS_USE_BOOTP */

#if (defined (RTEMS_SET_ETHERNET_ADDRESS))
  ethernet_address,               /* Ethernet hardware address */
#else
  NULL,                           /* Driver supplies hardware address */
#endif /* RTEMS_SET_ETHERNET_ADDRESS */
  0,				/* ignore broadcast */
  0,				/* mtu */
  0,				/* rbuf_count */
  0,				/* xbuf_count */
  0,				/* port */
#if defined(pc386)
  9,                      // irno -- isr vector
#elif defined(virtex)
  1,                      // irno -- isr vector
#else
  0,                      // irno -- isr vector
#endif
  0x81200000              // bpar -- hardware address (for virtex BSP)
};


/*
 * Network configuration
 */
struct rtems_bsdnet_config rtems_bsdnet_config = {
  &netdriver_config,
#if (defined (RTEMS_USE_BOOTP))
  rtems_bsdnet_do_bootp,
#else
  NULL,                /* do not use bootp */
#endif
  0,                   /* Default network task priority */
  0,                   /* Default mbuf capacity */
  0,                   /* Default mbuf cluster capacity */
  "rtems_ml405",       /* Host name */
  "extra.cea.fr",      /* Domain name */
#if defined(pc386)     /* QEMU network config */
  "10.0.2.2",          /* Gateway as seen from host PC when running qemu */
  "127.0.0.1",         /* Log host */
  {"10.0.2.3"},        /* Name servers(s) */
  {"10.0.2.3" },       /* NTP server(s) */
#elif (!defined (RTEMS_USE_BOOTP))
  "132.166.31.254",    /* Gateway */
  "127.0.0.1",		   /* Log host */
  {"132.167.198.4" },  /* Name server(s) */
  {"132.167.198.18" },    /* NTP server(s) */
#endif /* pc386 */
};

/*
 * For TFTP test application
 */
#if (defined (RTEMS_USE_BOOTP))
#define RTEMS_TFTP_TEST_HOST_NAME "BOOTP_HOST"
#define RTEMS_TFTP_TEST_FILE_NAME "BOOTP_FILE"
#else
#define RTEMS_TFTP_TEST_HOST_NAME "132.166.9.92"
#define RTEMS_TFTP_TEST_FILE_NAME "tftptest"
#endif

/*
 * For NFS test application
 * 
 * NFS server/path to mount and a directory to ls once mounted
 */
#define RTEMS_NFS_SERVER      "132.166.9.92"
#define RTEMS_NFS_SERVER_PATH "/home"
#define RTEMS_NFS_LS_PATH     "/media/nfstest"

#endif /* _RTEMS_NETWORKCONFIG_H_ */
