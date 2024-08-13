// =============================================================================
// Configs.cpp
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "CCfg/View/Configs.h"
#include "CCfg/Comment.h"

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture 1 definition
///////////////////////////////////////////////////////////////////////////////

Config1::Config1() : cfg("cfg")
{
	a_ = new Object("A", "*");
	a_x = new Scalar("x");
	a_y = new Scalar("y");
	a1 = new Object("A", "1");
	a1y = new Scalar("y");

	*a_x = -1.5;
	*a_y = 0;
	*a1y = 5;

	cfg.append(a_);
	cfg.append(a1);
	a_->append(a_x);
	a_->append(a_y);
	a1->append(a1y);
}

Config1::~Config1()
{

}

///////////////////////////////////////////////////////////////////////////////
// test fixture 2 definition
///////////////////////////////////////////////////////////////////////////////

Config2::Config2() : a("A")
{
	ab_ = new Object("B", "*");
	ab_c_ = new Object("C", "*");
	ab_c_d0 = new Object("D", "0");
	x = new Scalar("x");
	*x = 4.5;
	x->setPreComment(Comment("hello"));
	x->setPostComment(Comment("world"));

	ab0 = new Object("B", "0");
	ab0c0 = new Object("C", "0");
	ab0c0d_ = new Object("D", "*");
	y = new Scalar("y");
	*y = -2.3;

	ab0c0d0 = new Object("D", "0");
	z = new Scalar("z");
	*z = 0.1;

	ab0c0d0->setIncludeUrl("foo.xml");

	a.append(ab_);
	ab_->append(ab_c_);
	ab_c_->append(ab_c_d0);
	ab_c_d0->append(x);

	a.append(ab0);
	ab0->append(ab0c0);
	ab0c0->append(ab0c0d_);
	ab0c0d_->append(y);
	ab0c0->append(ab0c0d0);
	ab0c0d0->append(z);
}

Config2::~Config2()
{

}

///////////////////////////////////////////////////////////////////////////////
// test fixture 3 definition
///////////////////////////////////////////////////////////////////////////////

Config3::Config3() : cfg("cfg")
{
	a_ = new Object("A", "*");
	a_b0 = new Object("B", "0");
	a_b0c0 = new Object("C", "0");
	a_b0c0d0 = new Object("D", "0");
	a_b0c0d0x = new Scalar("x");
	a_b0c0d0y = new Scalar("y");
	a_b0c0d0z = new Scalar("z");
	*a_b0c0d0x = 4.5;
	*a_b0c0d0y = -0.5;
	*a_b0c0d0z = 0.5;

	a0 = new Object("A", "0");
	a0b_ = new Object("B", "*");
	a0b_c_ = new Object("C", "*");
	a0b_c_d0 = new Object("D", "0");
	a0b_c_d0y = new Scalar("y");
	a0b_c_d0w = new Scalar("w");
	*a0b_c_d0y = 2.5;
	*a0b_c_d0w = -6.5;

	a0b0 = new Object("B", "0");
	a0b0c0 = new Object("C", "0");
	a0b0c0d_ = new Object("D", "*");
	a0b0c0d_x = new Scalar("x");
	*a0b0c0d_x = 0;

	cfg.append(a_);
	a_->append(a_b0);
	a_b0->append(a_b0c0);
	a_b0c0->append(a_b0c0d0);
	a_b0c0d0->append(a_b0c0d0x);
	a_b0c0d0->append(a_b0c0d0y);
	a_b0c0d0->append(a_b0c0d0z);

	cfg.append(a0);
	a0->append(a0b_);
	a0b_->append(a0b_c_);
	a0b_c_->append(a0b_c_d0);
	a0b_c_d0->append(a0b_c_d0y);
	a0b_c_d0->append(a0b_c_d0w);

	a0->append(a0b0);
	a0b0->append(a0b0c0);
	a0b0c0->append(a0b0c0d_);
	a0b0c0d_->append(a0b0c0d_x);
}

Config3::~Config3()
{

}

///////////////////////////////////////////////////////////////////////////////
// test fixture 4 definition
///////////////////////////////////////////////////////////////////////////////

Config4::Config4() : cfg("cfg")
{
	a_ = new Object("A", "*");
	a_b_ = new Object("B", "*");
	a_b0 = new Object("B", "0");
	a_b0c_ = new Object("C", "*");
	a_b0c_u = new Scalar("u");


	a0 = new Object("A", "0");
	a0b_ = new Object("B", "*");
	a0b_c_ = new Object("C", "*");
	a0b_c0 = new Object("C", "0");
	a0b0 = new Object("B", "0");

	cfg.append(a_);
	a_->append(a_b_);
	a_->append(a_b0);
	a_b0->append(a_b0c_);
	a_b0c_->append(a_b0c_u);


	cfg.append(a0);
	a0->append(a0b_);
	a0b_->append(a0b_c_);
	a0b_->append(a0b_c0);
	a0->append(a0b0);
}

Config4::~Config4()
{

}

///////////////////////////////////////////////////////////////////////////////
// test fixture 5 definition
///////////////////////////////////////////////////////////////////////////////

Config5::Config5() : cfg("cfg")
{
	a_ = new Object("A", "*");
	a_b0 = new Scalar("b", "0");
	a_b1 = new Scalar("b", "1");
	a_c0 = new Scalar("c", "0");
	a_d0 = new Scalar("d", "0");
	a0 = new Object("A", "0");
	a0b2 = new Scalar("b", "2");
	a0c1 = new Scalar("c", "1");
	a0c2 = new Scalar("c", "2");
	a0d1 = new Scalar("d", "1");

	cfg.append(a_);
	a_->append(a_b0);
	a_->append(a_b1);
	a_->append(a_c0);
	a_->append(a_d0);
	cfg.append(a0);
	a0->append(a0b2);
	a0->append(a0c1);
	a0->append(a0c2);
	a0->append(a0d1);
}

Config5::~Config5()
{

}

}
