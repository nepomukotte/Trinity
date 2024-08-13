/*
 * AttributeClass.cpp
 *
 *  Created on: 26 juin 2009
 *      Author: kbenzi
 */

#include "AttributeClass.h"
#include "CCfg/Attribute.h"
#include "CCfg/Range.h"
#include "CCfg/Interval.h"
#include "CCfg/Enum.h"
#include "CCfg/Path.h"
#include <memory>

namespace CCfg
{

AttributeClass::AttributeClass()
{

}

/**
 * Attributeclass constructor from an Attribute object
 * @param attr Attribute for building the Attribute Class
 */
AttributeClass::AttributeClass(const Attribute* attr)
:	path(attr->getShortClassPath().toString()),
	name(attr->getName()),
	unit(attr->isObject() ? "" : attr->getUnit()),
	attrTypeId(attr->getType())
{
	Range* range = attr->isObject() ? 0 : attr->getRange();
	if (range != 0)
	{
		// Finds the type (int, real, float, hex, string)
		ValueTypeId valueType = attr->getValueTypeId();

		Interval* interv = dynamic_cast<Interval*>(range);
		if (interv != 0) // Range is an Interval
		{
			switch(valueType)
			{
			case Id<Int>::val:
				valueRange.reset(new ValueRangeInterval<Int>(interv));
				break;
			case Id<Real>::val:
				valueRange.reset(new ValueRangeInterval<Real>(interv));
				break;
			case Id<Hex>::val:
				valueRange.reset(new ValueRangeInterval<Hex>(interv));
				break;
			case Id<Bool>::val:
				valueRange.reset(new ValueRangeInterval<Bool>(interv));
				break;
			case Id<String>::val:
				valueRange.reset(new ValueRangeInterval<String>(interv));
				break;
			default:
				break;
			}
		}
		Enum* enumRange = dynamic_cast<Enum*>(range);
		if (enumRange != 0) // Range is an Enum
		{
			//Find the enum type
			ValueTypeId valueType = enumRange->getValues().begin()->getType();
			// Find the type (int, real, float, hex, string)
			switch(valueType)
			{
			case Id<Int>::val:
				valueRange.reset(new ValueRangeEnum<Int>(enumRange));
				break;
			case Id<Real>::val:
				valueRange.reset(new ValueRangeEnum<Real>(enumRange));
				break;
			case Id<Hex>::val:
				valueRange.reset(new ValueRangeEnum<Hex>(enumRange));
				break;
			case Id<Bool>::val:
				valueRange.reset(new ValueRangeEnum<Bool>(enumRange));
				break;
			case Id<String>::val:
				valueRange.reset(new ValueRangeEnum<String>(enumRange));
				break;
			default:
				break;
			}
		}
	}
}


}//namespace CCfg


