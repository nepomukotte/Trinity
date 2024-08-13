#ifndef __SMB_SLAVE_H__
#define __SMB_SLAVE_H__

/*****************************************************************************
* Configuration parameters.
* Edit the following lines to customize the SMBus driver.
*****************************************************************************/
/*
 * Uncomment one of the following two lines to support PEC.
 * If both lines are commented out, no PEC support is included.
 */
//#define SMB_USE_PEC_LOOKUP          //!< Use CRC lookup table for PEC
#define SMB_USE_PEC_CALCULATION   //!< Use CRC calculation for PEC

//! The 7 bit slave address of this device.
#define SMB_OWN_ADDRESS       0x10 // Reserved address. Not for production.

/*!
 *  Maximum number of data bytes received for Block write and Block write,
 *  block read process call. (Max value is 32).
 */
#define SMB_RX_MAX_LENGTH       32

/*
 *  Maximum number of data bytes transmitted for Block read and Block write,
 *  block read process call. (Max value is 32).
 */
#define SMB_TX_MAX_LENGTH       32

/*****************************************************************************
* End of configuration parameter section
*****************************************************************************/


#ifdef SMB_USE_PEC_LOOKUP
#define SMB_SUPPORT_PEC
#endif

#ifdef SMB_USE_PEC_CALCULATION
#define SMB_SUPPORT_PEC
#endif

//! Length of command code.
#define SMB_COMMAND_CODE_LENGTH   1

//! Length of byte count.
#define SMB_BYTE_COUNT_LENGTH     1

//! Length of PEC.
#define SMB_PEC_LENGTH            1

//! Length of receive buffer, must be large enough to include control bytes.
#ifdef SMB_SUPPORT_PEC
#define SMB_RX_BUFFER_LENGTH    (SMB_COMMAND_CODE_LENGTH + SMB_BYTE_COUNT_LENGTH + SMB_RX_MAX_LENGTH + SMB_PEC_LENGTH)
#else
#define SMB_RX_BUFFER_LENGTH    (SMB_COMMAND_CODE_LENGTH + SMB_BYTE_COUNT_LENGTH + SMB_RX_MAX_LENGTH)
#endif

/*!
 * Length of transmit buffer, must be large enough to include control bytes.
 * Room for PEC is not needed, since it will never be placed in the transmit buffer.
 */
#define SMB_TX_BUFFER_LENGTH    (SMB_BYTE_COUNT_LENGTH + SMB_TX_MAX_LENGTH)

//! Value of write bit appended after slave address in SMBus communication.
#define SMB_WRITE                       0

//! Value of read bit appended after slave address in SMBus communication.
#define SMB_READ                        1

//! Value of slave address with write bit appended (used for PEC calculation/lookup).
#define SMB_OWN_ADDRESS_W               ((SMB_OWN_ADDRESS << 1) | SMB_WRITE)

//! Value of slave address with reaad bit appended (used for PEC calculation/lookup).
#define SMB_OWN_ADDRESS_R               ((SMB_OWN_ADDRESS << 1) | SMB_READ)

//! The CRC polynome used in PEC calculation.
#define SMB_CRC_POLYNOME                0x07

/*! \brief Macro that checks the length and PEC of a message.
 *
 *  This macro can be used to check if a message is of length n. If it is
 *  of length n+1 and PEC is enabled, the PEC is checked.
 *
 *  Returns 0 (FALSE) if message is of wrong length or on PEC error.
 */
#ifdef SMB_SUPPORT_PEC
#define RX_COUNT_AND_PEC_CHECK(n)  (smb->rxCount == n || (smb->rxCount == (n + 1) && smb->pec == 0))
#else
#define RX_COUNT_AND_PEC_CHECK(n)  (smb->rxCount == n)
#endif




#define SMB_STATE_IDLE                  0x00    //!< Idle state flag.
#define SMB_STATE_READ_REQUESTED        0x01    //!< Read requested flag.
#define SMB_STATE_WRITE_REQUESTED       0x02    //!< Write requested flag.
#define SMB_STATE_WRITE_READ_REQUESTED  0x03    //!< Write, then read requested flag.


#define TRUE    1
#define FALSE   0


// Function prototypes.
void SMBusInit(void);
void SMBEnable(void);
void SMBDisable(void);
unsigned char SMBError(void);
void SMBcrc(unsigned char data);


/*! \brief Collection of all variables used by SMBus driver.
 *
 *  The SMBData struct contains all the variables used internally by the SMBus slave.
 */
typedef struct SMBData
{
    unsigned char txLength;                         //!< Transmit length.
    unsigned char txCount;                          //!< Transmit counter.
    unsigned char rxCount;                          //!< Receive counter.
    unsigned char state:2;                          //!< SMBus driver state flag.
    unsigned char volatile enable : 1;              //!< Enable ACK on requests.
    unsigned char volatile error : 1;               //!< Error flag.
#ifdef SMB_SUPPORT_PEC
    unsigned char pec;                              //!< PEC of message in progress.
#endif //SMB_SUPPORT_PEC
    unsigned char rxBuffer[SMB_RX_BUFFER_LENGTH];   //!< Receive buffer.
    unsigned char txBuffer[SMB_TX_BUFFER_LENGTH];   //!< Transmit buffer.
} SMBData;


/*! \brief Struct with global flags related to the SMBus driver.
 *
 *  The SMBGlobalFlags struct contains global flags related to the SMBus driver.
 */
/*
typedef struct SMBGlobalFlags
{
    unsigned char enable : 1;                       //!< Enable ACK on requests
    unsigned char error : 1;                        //!< Error flag
} SMBGlobalFlags;
*/

#endif
