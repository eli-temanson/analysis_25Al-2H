#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include <math.h>
#include "catima/catima.h"
#include "catima/reactions.h"
#include "testutils.h"
using namespace std;
using catima::reaction_rate;
using catima::nonreaction_rate;
    TEST_CASE("reaction"){
        catima::Projectile proj{12,6,6,870};
        catima::Projectile proj2{238,92,92,500};
        auto c = catima::get_material(6);
        auto h = catima::get_material(1);
        catima::Material water({{0,8,2},{0,1,1}},1.0);
        c.thickness(2.0);
        double r,r2;

        r= catima::nonreaction_rate(proj, c);
        CHECK(r == approx(0.92,0.02));

        catima::Layers l;
        l.add(c);
        l.add(c);
        l.add(c);

        auto res = catima::calculate(proj,l);
        CHECK(res.total_result.sp == approx(r*r*r,0.02));

        c.thickness(6.0);
        r= catima::nonreaction_rate(proj, c);
        CHECK(res.total_result.sp == approx(r,0.001));
        
        c.thickness(0.0);
        r= catima::nonreaction_rate(proj, c);
        CHECK(r == 1.0);
        proj.T = 0.0;
        c.thickness(6);
        r= catima::nonreaction_rate(proj, c);
        CHECK(r == -1.0);

        proj.T=870;
        water.thickness(1);
        r= catima::nonreaction_rate(proj2, water);
        r2= catima::nonreaction_rate(proj, water);
        CHECK( (r > 0 && r<1.0) );
        CHECK( (r2 > 0 && r2<1.0) );
        CHECK( r2>r );
    }
    TEST_CASE("production"){
        catima::Projectile proj{12,6,6,870};
        auto c = catima::get_material(6);
        c.thickness(0.1);
        double r,r2;
        
        double cs = 70;
        double rcsi = 870;
        double rcso = 850;
        
        r = 0.0001*cs*c.number_density_cm2();
        r2 = catima::production_rate(cs,rcsi,rcso, c);
        CHECK(r==approx(r2).R(0.01));
        
        r2 = catima::production_rate(cs,2,1, c);
        CHECK(r==approx(r2).R(0.001));
        
        r = catima::production_rate(cs,870,870, c);
        r2 = catima::production_rate(cs,870,860, c);
        CHECK(r==approx(r2).R(0.001));
        c.thickness(2.0);
        r = catima::production_rate(cs,870,870, c);
        r2 = catima::production_rate(cs,870,860, c);
        CHECK(r==approx(r2).R(0.001));
    }
