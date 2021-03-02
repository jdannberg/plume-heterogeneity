/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/initial_composition/interface.h>
#include <aspect/boundary_composition/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class ChemicalHeterogeneity : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the compositional fields.
         */
        std::unique_ptr<Functions::ParsedFunction<dim> > function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        unsigned int n_blobs;
        double minimum_blob_radius;
        double maximum_blob_radius;
    };

  }

  namespace BoundaryComposition
  {
    using namespace dealii;

    /**
     * A class that implements boundary composition based on a functional
     * description provided in the input file.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class ChemicalHeterogeneity : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the boundary composition as a function of position and time.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the compositional fields.
         */
        std::unique_ptr<Functions::ParsedFunction<dim> > function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        std::vector<Point<dim> > blob_centers;
        std::vector<double> blob_radii;
        double end_time;
        double blob_spacing;
        double minimum_blob_radius;
        double maximum_blob_radius;
        double average_velocity;
    };
  }
}

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    double
	ChemicalHeterogeneity<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // this initial condition only makes sense if the geometry is a
      // Box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != nullptr,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(1);
      Point<dim> blob_center;
      double blob_radius;
      double blob_present = 0;
      const Point<dim> extents = geometry->get_extents();

      for (unsigned int n=0; n<n_blobs; ++n)
        {
    	  do
    	  {
            for (unsigned int d=0; d<dim; ++d)
              blob_center[d] = (double)(std::rand() % int(extents[d]));
    	  } while(!geometry->point_is_in_domain(blob_center) || function->value(blob_center,n_comp) == 0.0);

          blob_radius = minimum_blob_radius + std::rand() % int(maximum_blob_radius - minimum_blob_radius);

          if (position.distance(blob_center) < blob_radius)
        	blob_present = 1.0;
        }

      return blob_present;
    }


    template <int dim>
    void
	ChemicalHeterogeneity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Chemical Heterogeneity");
        {
          /**
           * Choose the coordinates to evaluate the maximum refinement level
           * function. The function can be declared in dependence of depth,
           * cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2D/3D "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);

          prm.declare_entry ("Number of blobs", "10",
                             Patterns::Double (0),
                             "Number of spherical blobs with a different composition.");
          prm.declare_entry ("Minimum blob radius", "10000",
                             Patterns::Double (0),
                             "Minimum radius of spherical blobs with a different composition. Units: $\\si{m}$.");
          prm.declare_entry ("Maximum blob radius", "25000",
                             Patterns::Double (0),
                             "Maximum radius of spherical blobs with a different composition. Units: $\\si{m}$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
	ChemicalHeterogeneity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Chemical Heterogeneity");
        {
          coordinate_system   = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
          n_blobs             = prm.get_double ("Number of blobs");
          minimum_blob_radius = prm.get_double ("Minimum blob radius");
          maximum_blob_radius = prm.get_double ("Maximum blob radius");
        }

        try
          {
            function
              = std_cxx14::make_unique<Functions::ParsedFunction<dim>>(this->n_compositional_fields());
            function->parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Initial composition model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'\n"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }



  namespace BoundaryComposition
  {

    template <int dim>
    double
	ChemicalHeterogeneity<dim>::
    boundary_composition (const types::boundary_id /*boundary_indicator*/,
                          const Point<dim> &position,
                          const unsigned int /*compositional_field*/) const
    {
      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(1);

      Point<dim> blob_position(position);
      blob_position(dim-1) += average_velocity * this->get_time();

      for (unsigned int n=0; n<blob_centers.size(); ++n)
        {
          if (blob_position.distance(blob_centers[n]) < blob_radii[n])
        	return 1.0;
        }

      return 0.0;
    }


    template <int dim>
    void
	ChemicalHeterogeneity<dim>::initialize ()
    {
      // this boundary condition only makes sense if the geometry is a
      // Box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != nullptr,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      // find out how many blobs we will need
      const unsigned int n_blobs = average_velocity / blob_spacing * end_time;
      blob_centers.resize(n_blobs);
      blob_radii.resize(n_blobs);

      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(1);
      const Point<dim> extents = geometry->get_extents();

      for (unsigned int n=0; n<n_blobs; ++n)
        {
          blob_centers[n](0) = (double)(std::rand() % int(extents[0]));
          blob_centers[n](1) = (double)(std::rand() % int(average_velocity * end_time));
          blob_radii[n] = minimum_blob_radius + std::rand() % int(maximum_blob_radius - minimum_blob_radius);
        }

      return;
    }



    template <int dim>
    void
	ChemicalHeterogeneity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Chemical Heterogeneity");
        {
          prm.declare_entry ("Blob spacing", "15000",
                             Patterns::Double (0),
                             "Average vertical spacing between the blobs that flow in as a boundary "
                             "composition. Location of the blobs is random.");
          prm.declare_entry ("Average velocity", "0.1",
                             Patterns::Double (0),
                             "Average upwelling velocity, which is used to compute the spacing of "
                             "the inflowing blobs.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
	ChemicalHeterogeneity<dim>::parse_parameters (ParameterHandler &prm)
    {
      // read end time from parameter file. if it is to be interpreted
      // in years rather than seconds, then do the conversion
      end_time = prm.get_double ("End time");
      if (this->convert_output_to_years())
        end_time *= year_in_seconds;

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Chemical Heterogeneity");
        {
          minimum_blob_radius = prm.get_double ("Minimum blob radius");
          maximum_blob_radius = prm.get_double ("Maximum blob radius");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Chemical Heterogeneity");
        {
          blob_spacing = prm.get_double ("Blob spacing");

          average_velocity = prm.get_double ("Average velocity");
          if (this->convert_output_to_years() == true)
        	  average_velocity /= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(ChemicalHeterogeneity,
                                              "chemical heterogeneity",
                                              ".")
  }

  namespace BoundaryComposition
  {
  ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(ChemicalHeterogeneity,
                                             "chemical heterogeneity",
                                             ".")
  }
}
