/* =====================================================================================================================
 * DaqControl.ice
 * ---------------------------------------------------------------------------------------------------------------------
 * mdaq::DaqControl remote interface definition.
 * Created on: 4 May 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mdaq_DaqControl_ice_INCLUDED
#define mdaq_DaqControl_ice_INCLUDED

#include "mdaq/utl/CmdException.ice"

module mdaq
{
	interface DaqCtrlNode
	{
		/**
		  * Creates a sender for data from the circular buffer.
		  * @param dataRouterType Protocol to use for sending data (Ice, Tcp, ZBuf or Fdt).
		  * @param dataRouterName Endpoint of the data receiver.
		  * @throws ::mdaq::utl::CmdException Raised in failure to create the data sender.
		 **/
	    void connect(string dataRouterType, string dataRouterName) throws ::mdaq::utl::CmdException;

		/** Closes any open connection. **/
	    idempotent void disconnect();
	    /**
	      * Sets data sender policy.
	      * @param enable Whether to sent all data immediately or to optimize data chunk size.
	      * @throws ::mdaq::utl::CmdException Raised if no data sender currently exists.
	      **/
	    void setAlwaysFlushData(bool enable) throws ::mdaq::utl::CmdException;
	    
		/**
	      * Selects which of the circular buffers should be read.
	      * @param mask Mask of enabled buffers (1 bit per buffer).
	      **/
	    void setCircularBuffersEnabled(int mask);

		/** Establishes connection from data sender to data receiver and enables interrupts. **/
	    void daqStart() throws ::mdaq::utl::CmdException;
	    
	    /** Disables interrupts and closes data link. **/
	    void daqStop() throws ::mdaq::utl::CmdException;

		void shutdown();
	};  
};

#endif /* mdaq_DaqControl_ice_INCLUDED */
