/* =============================================================================
 * ServerDevices.cpp
 * -----------------------------------------------------------------------------
 * Created on: 27 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the software for GET.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#include "mdaq/hw/server/DynDevice.hpp"
#include "mdaq/hw/server/DeviceFactory.h"
#include "mdaq/hw/NestedStoragePolicy.hpp"
#include "mdaq/hw/DeviceStoragePolicy.hpp"
#include "mdaq/hw/ControllerDeviceSimulator.hpp"

using namespace mdaq::hw;
using namespace mdaq::hw::server;

REGISTER_DEVICE_TYPE(uint32_t, NestedStoragePolicy, NonControllerDevice, Nested)
REGISTER_DEVICE_TYPE(uint32_t, NestedStoragePolicy, NonControllerDevice, MemBus)
REGISTER_DEVICE_TYPE(uint64_t, NestedStoragePolicy, NonControllerDevice, AGetBus)
REGISTER_DEVICE_TYPE(uint64_t, NestedStoragePolicy, NonControllerDevice, AsAdBus)
REGISTER_DEVICE_TYPE(uint32_t, DeviceStoragePolicy, ControllerDeviceSimulator, LocalMap)

