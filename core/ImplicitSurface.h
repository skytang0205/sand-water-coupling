#pragma once

#include "Surface.h"

namespace Pivot {

	class ImplicitSphere : public Surface {
	public:

		ImplicitSphere(Vector3d const &center, double radius) : m_Center { center }, m_Radius { radius } { }

		virtual Vector3d ClosestPositionOf(Vector3d const &pos) const override { return m_Center + ClosestNormalOf(pos) * m_Radius; }
		virtual Vector3d ClosestNormalOf  (Vector3d const &pos) const override { return (pos - m_Center).normalized(); }
		virtual double   SignedDistanceTo (Vector3d const &pos) const override { return (pos - m_Center).norm() - m_Radius; }

		virtual std::pair<Vector3d, Vector3d> GetCornersOfAABB() const override {
			return { m_Center - Vector3d::Ones() * m_Radius, m_Center + Vector3d::Ones() * m_Radius };
		}

	private:
		Vector3d m_Center;
		double   m_Radius;
	};

	class ImplicitBox : public Surface {
	public:

		ImplicitBox(Vector3d const &minCorner, Vector3d const &lengths) : m_Center { minCorner + lengths / 2 }, m_HalfLengths { lengths / 2 } { }

		virtual Vector3d ClosestNormalOf(Vector3d const &pos) const override {
			Vector3d const phi = (pos - m_Center).cwiseAbs() - m_HalfLengths;
			Vector3d normal;
			if ((phi.array() <= 0).all()) {
				int axis;
				phi.maxCoeff(&axis);
				normal = Vector3d::Unit(axis);
			} else {
				normal = phi.cwiseMax(0);
			}
			return normal.cwiseProduct((pos - m_Center).cwiseSign()).normalized();
		}

		virtual double SignedDistanceTo(Vector3d const &pos) const override {
			Vector3d const phi = (pos - m_Center).cwiseAbs() - m_HalfLengths;
			if ((phi.array() <= 0).all()) {
				return phi.maxCoeff();
			} else {
				return phi.cwiseMax(0).norm();
			}
		}

		virtual std::pair<Vector3d, Vector3d> GetCornersOfAABB() const override {
			return { m_Center - m_HalfLengths, m_Center + m_HalfLengths };
		}

	private:
		Vector3d m_Center;
		Vector3d m_HalfLengths;
	};

	class ImplicitPlane : public Surface {
	public:
		ImplicitPlane(Vector3d const &position, Vector3d const &direction) : m_Position(position), m_Normal(direction.normalized()) { }

		virtual Vector3d ClosestNormalOf (Vector3d const &pos) const override { return m_Normal; }
		virtual double   SignedDistanceTo(Vector3d const &pos) const override { return (pos - m_Position).dot(m_Normal); }

		virtual std::pair<Vector3d, Vector3d> GetCornersOfAABB() const override {
			return { Vector3d::Zero(), Vector3d::Zero() };
		}

	private:
		Vector3d const m_Position;
		Vector3d const m_Normal;
	};

	class ImplicitEllipsoid : public Surface {
	public:

		ImplicitEllipsoid(Vector3d const &center, Vector3d const &semiAxels) : m_Center(center), m_SemiAxels(semiAxels) { }

		virtual Vector3d ClosestPositionOf(Vector3d const &pos) const override { return pos / pos.cwiseQuotient(m_SemiAxels).norm(); } // not accurate solution
		virtual Vector3d ClosestNormalOf  (Vector3d const &pos) const override { return (pos - ClosestPositionOf(pos)).normalized() * (Surrounds(pos) ? -1 : 1); }
		virtual double   DistanceTo       (Vector3d const &pos) const override { return (pos - ClosestPositionOf(pos)).norm(); }
		virtual double   SignedDistanceTo (Vector3d const &pos) const override { return DistanceTo(pos) * (Surrounds(pos) ? -1 : 1); }
		virtual bool     Surrounds        (Vector3d const &pos) const override { return pos.cwiseQuotient(m_SemiAxels).squaredNorm() <= 1; }

		virtual std::pair<Vector3d, Vector3d> GetCornersOfAABB() const override {
			return { m_Center - m_SemiAxels, m_Center + m_SemiAxels };
		}
		
	private:
		Vector3d m_Center;
		Vector3d m_SemiAxels;
	};
}
