#include <iostream>
#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/Geometry>
#include <CImg.h>
#include <vector>
#include <string>

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

//typedef std::complex<double> complexd;
const double pi = 3.141562;

namespace Drehkristall
{

    class Atom
    {
        public:
        Atom(const double r, const double z, const Vector3d &vek)
        : radius(r), position(vek), elektronenzahl(z)
        {
            radiusquadrat=radius*radius;
        }

        Vector3d position;
        double radius;
        double elektronenzahl;
        double radiusquadrat;

        double Atomformfaktor(const Vector3d &dk)
        {
            return(1.0);//elektronenzahl/(1.0+radiusquadrat*k.cross(k0).dot(k.cross(k0))));
        }
    };

    class BasisKubisch
    {
        public:
        BasisKubisch(const double gitterkonstante)
        {
            a=gitterkonstante;
        }

        void AddAtom(Atom &atom)
        {
            atome.push_back(atom);
        }

        double getPhase(const Vector3d &dk)
        {
            double phase=0;
            std::vector<Atom>::iterator it;
            for(it=atome.begin();it!=atome.end();it++)
            {
                phase+=it->Atomformfaktor(dk)*cos(dk.dot(it->position)*a);
            }
            return(phase);
        }
        Vector3d getGitterVektor(unsigned int num)
        {

        }

        Vector3d getReziprokenGitterVektor(unsigned int num)
        {
            switch(num)
            {
                case 0:
                return(Vector3d(2.0*pi/a,0,0));
                case 1:
                return(Vector3d(0,2.0*pi/a,0));
                case 2:
                return(Vector3d(0,0,2.0*pi/a));
            }
        }

        double a;
        std::vector<Drehkristall::Atom> atome;
    };

    class dkfilm
    {
        public:
        dkfilm( const unsigned int w, const unsigned int h, const double r, const double heightmm, const double mi=1.0)
        :   width(w), height(h), maxintens(mi), film(w,h,1,4,0), radius(r),
            pixelsperradian(w/(2.0*pi)), pixelspermm(h/heightmm)
        {
            white[0] = white[1] = white[2] = white[3] = 255;
        }

        ~dkfilm()
        {

        }

        void reflex( double alpha, double beta, double intensity, bool middlebeam=true)
        {
            if(!middlebeam)
                alpha+=pi;
            if(intensity > maxintens)
                intensity = maxintens;
            const unsigned char intenscolor[] = {   static_cast<unsigned char>(255*intensity/maxintens),
                                                    static_cast<unsigned char>(255*intensity/maxintens),
                                                    static_cast<unsigned char>(255*intensity/maxintens), 255 };
            while(alpha<0)
                alpha+=2.0*pi;
            while(alpha>2.0*pi)
                alpha-=2.0*pi;
            while(beta<0)
                beta+=2.0*pi;
            while(beta>2.0*pi)
                beta-=2.0*pi;
            film.draw_circle(alpha*pixelsperradian, height/2-tan(beta)*radius, 2, intenscolor, 0, 0.5);
        }

        void save(const char *str)
        {
            film.save(str);
        }

        private:

        const double width,height;
        const double maxintens;
        cimg_library::CImg<unsigned char> film;
        const double radius;
        const double pixelsperradian;
        const double pixelspermm;

        unsigned char white[4];
    };
}
int main()
{
    const unsigned int width = 1080;
    const unsigned int height = 400;
    const double radsperpixel = 0.005817764;
    const double mmperpixel = 0.25;

    Drehkristall::BasisKubisch blub(0.5604);
    Drehkristall::Atom Na(0.098,11,Vector3d(0,0,0));
    Drehkristall::Atom Cl(0.181,17,Vector3d(0.5,0,0));
    blub.AddAtom(Na);
    blub.AddAtom(Cl);

    const double lambda=0.154178;
    const double k0_len=2*pi/lambda;
    const Vector3d k0_ur(0,0,1);
    const Vector3d k0up_ur(0,1,0);
    const Vector3d k0left_ur=k0_ur.cross(k0up_ur);

    //Gitter gitter(1,1,1,&blub,1);

    const double radius=287500000.0;
    const Vector3d rotationsachse(0,1,0);

    const signed int numh=40;
    const signed int numk=40;
    const signed int numl=40;

    Vector3d rgvektor0 = blub.getReziprokenGitterVektor(0);
    Vector3d rgvektor1 = blub.getReziprokenGitterVektor(1);
    Vector3d rgvektor2 = blub.getReziprokenGitterVektor(2);

    Drehkristall::dkfilm film(1080,400,300.0,40.0,1.0);

    int x,y,z;

    std::cout << -numh/2 << std::endl;
    for(x=-numh/2;x<numh/2;x++)
    {
        for(y=-numk/2;y<numk/2;y++)
        {
            for(z=-numl/2;z<numl/2;z++)
            {
                Vector3d k = x*rgvektor0+y*rgvektor1+z*rgvektor2;
                double klen = k.norm();
                if(klen>2.0*k0_len)
                    continue;

                double alpha = acos(klen/(2.0*k0_len));

                double beta = atan(k.z()/k0_len);

                double intens = blub.getPhase(k);

                std::cout << x << "," << y << "," << z << ":" << intens << std::endl;
                std::cout << alpha << "," << beta << std::endl;

                film.reflex(2.0*alpha, beta, intens);
                film.reflex(-2.0*alpha, beta, intens);
            }
        }
    }

    film.save("test.png");

    return 0;
}
