using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeoAPI.CoordinateSystems;
using GeoAPI.Geometries;
using GeoAPI.CoordinateSystems.Transformations;

namespace ProjNet
{
    public class Gcj02toWgs84:GeoAPI.CoordinateSystems.Transformations.ICoordinateTransformation
    {
        public Gcj02toWgs84()
            : base()
        {
            _SourceCS = ProjNet.CoordinateSystems.GeocentricCoordinateSystem.WGS84;
            _TargetCS = ProjNet.CoordinateSystems.ProjectedCoordinateSystem.WebMercator;
            _TransformType = GeoAPI.CoordinateSystems.Transformations.TransformType.Other;
            _Remarks = "HT zhao & JG li";
            _Name = "Gcj02toWgs84";
            _MathTransform = null;
        }
        #region ICoordinateTransformation Members

        private string _AreaOfUse;
        /// <summary>
        /// Human readable description of domain in source coordinate system.
        /// </summary>		
        public string AreaOfUse
        {
            get { return _AreaOfUse; }
        }

        private string _Authority;
        /// <summary>
        /// Authority which defined transformation and parameter values.
        /// </summary>
        /// <remarks>
        /// An Authority is an organization that maintains definitions of Authority Codes. For example the European Petroleum Survey Group (EPSG) maintains a database of coordinate systems, and other spatial referencing objects, where each object has a code number ID. For example, the EPSG code for a WGS84 Lat/Lon coordinate system is ?326?
        /// </remarks>
        public string Authority
        {
            get { return _Authority; }
        }

        private long _AuthorityCode;
        /// <summary>
        /// Code used by authority to identify transformation. An empty string is used for no code.
        /// </summary>
        /// <remarks>The AuthorityCode is a compact string defined by an Authority to reference a particular spatial reference object. For example, the European Survey Group (EPSG) authority uses 32 bit integers to reference coordinate systems, so all their code strings will consist of a few digits. The EPSG code for WGS84 Lat/Lon is ?326?</remarks>
        public long AuthorityCode
        {
            get { return _AuthorityCode; }
        }

        private IMathTransform _MathTransform;
        /// <summary>
        /// Gets math transform.
        /// </summary>
        public IMathTransform MathTransform
        {
            get { return _MathTransform; }
        }

        private string _Name;
        /// <summary>
        /// Name of transformation.
        /// </summary>
        public string Name
        {
            get { return _Name; }
        }

        private string _Remarks;
        /// <summary>
        /// Gets the provider-supplied remarks.
        /// </summary>
        public string Remarks
        {
            get { return _Remarks; }
        }

        private ICoordinateSystem _SourceCS;
        /// <summary>
        /// Source coordinate system.
        /// </summary>
        public ICoordinateSystem SourceCS
        {
            get { return _SourceCS; }
        }

        private ICoordinateSystem _TargetCS;
        /// <summary>
        /// Target coordinate system.
        /// </summary>
        public ICoordinateSystem TargetCS
        {
            get { return _TargetCS; }
        }

        private TransformType _TransformType;
        /// <summary>
        /// Semantic type of transform. For example, a datum transformation or a coordinate conversion.
        /// </summary>
        public TransformType TransformType
        {
            get { return _TransformType; }
        }

        #endregion
    }

    public class Wgs84toGcj02 : GeoAPI.CoordinateSystems.Transformations.ICoordinateTransformation
    {
        public Wgs84toGcj02()
            : base()
        {
            _SourceCS = ProjNet.CoordinateSystems.GeocentricCoordinateSystem.WGS84;
            _TargetCS = ProjNet.CoordinateSystems.ProjectedCoordinateSystem.WebMercator;
            _TransformType = GeoAPI.CoordinateSystems.Transformations.TransformType.Other;
            _Remarks = "HT zhao & JG li";
            _Name = "Gcj02toWgs84";
            _MathTransform = null;
        }


        #region ICoordinateTransformation Members

        private string _AreaOfUse;
        /// <summary>
        /// Human readable description of domain in source coordinate system.
        /// </summary>		
        public string AreaOfUse
        {
            get { return _AreaOfUse; }
        }

        private string _Authority;
        /// <summary>
        /// Authority which defined transformation and parameter values.
        /// </summary>
        /// <remarks>
        /// An Authority is an organization that maintains definitions of Authority Codes. For example the European Petroleum Survey Group (EPSG) maintains a database of coordinate systems, and other spatial referencing objects, where each object has a code number ID. For example, the EPSG code for a WGS84 Lat/Lon coordinate system is ?326?
        /// </remarks>
        public string Authority
        {
            get { return _Authority; }
        }

        private long _AuthorityCode;
        /// <summary>
        /// Code used by authority to identify transformation. An empty string is used for no code.
        /// </summary>
        /// <remarks>The AuthorityCode is a compact string defined by an Authority to reference a particular spatial reference object. For example, the European Survey Group (EPSG) authority uses 32 bit integers to reference coordinate systems, so all their code strings will consist of a few digits. The EPSG code for WGS84 Lat/Lon is ?326?</remarks>
        public long AuthorityCode
        {
            get { return _AuthorityCode; }
        }

        private IMathTransform _MathTransform;
        /// <summary>
        /// Gets math transform.
        /// </summary>
        public IMathTransform MathTransform
        {
            get { return _MathTransform; }
        }

        private string _Name;
        /// <summary>
        /// Name of transformation.
        /// </summary>
        public string Name
        {
            get { return _Name; }
        }

        private string _Remarks;
        /// <summary>
        /// Gets the provider-supplied remarks.
        /// </summary>
        public string Remarks
        {
            get { return _Remarks; }
        }

        private ICoordinateSystem _SourceCS;
        /// <summary>
        /// Source coordinate system.
        /// </summary>
        public ICoordinateSystem SourceCS
        {
            get { return _SourceCS; }
        }

        private ICoordinateSystem _TargetCS;
        /// <summary>
        /// Target coordinate system.
        /// </summary>
        public ICoordinateSystem TargetCS
        {
            get { return _TargetCS; }
        }

        private TransformType _TransformType;
        /// <summary>
        /// Semantic type of transform. For example, a datum transformation or a coordinate conversion.
        /// </summary>
        public TransformType TransformType
        {
            get { return _TransformType; }
        }

        #endregion
    }

    public class Gcj02toWgs84Math : CoordinateSystems.Projections.MapProjection
    {
        public Gcj02toWgs84Math():base(null)
        { 

        }

        private double transformLat(double x, double y)
        {
            double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * Math.Sqrt(Math.Abs(x));
            ret += (20.0 * Math.Sin(6.0 * x * Math.PI) + 20.0 * Math.Sin(2.0 * x * Math.PI)) * 2.0 / 3.0;
            ret += (20.0 * Math.Sin(y * Math.PI) + 40.0 * Math.Sin(y / 3.0 * Math.PI)) * 2.0 / 3.0;
            ret += (160.0 * Math.Sin(y / 12.0 * Math.PI) + 320 * Math.Sin(y * Math.PI / 30.0)) * 2.0 / 3.0;
            return ret;
        }

        private double transformLon(double x, double y)
        {
            double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * Math.Sqrt(Math.Abs(x));
            ret += (20.0 * Math.Sin(6.0 * x * Math.PI) + 20.0 * Math.Sin(2.0 * x * Math.PI)) * 2.0 / 3.0;
            ret += (20.0 * Math.Sin(x * Math.PI) + 40.0 * Math.Sin(x / 3.0 * Math.PI)) * 2.0 / 3.0;
            ret += (150.0 * Math.Sin(x / 12.0 * Math.PI) + 300.0 * Math.Sin(x / 30.0 * Math.PI)) * 2.0 / 3.0;
            return ret;
        }
			
		/// <summary>
		/// Converts coordinates in decimal degrees to projected meters.
		/// </summary>
		/// <param name="lonlat">The point in decimal degrees.</param>
		/// <returns>Point in projected meters</returns>
        protected override double[] RadiansToMeters(double[] lonlat)
		{
            var lon = lonlat[0];		   
            var lat = lonlat[1];    

            double ee = 0.00669342162296594323;
            double a = 6378245.0;
            double dLat = transformLat(lon - 105.0, lat - 35.0);
            double dLon = transformLon(lon - 105.0, lat - 35.0);
            double radLat = lat / 180.0 * Math.PI;
            double magic = Math.Sin(radLat);
            magic = 1 - ee * magic * magic;
            double sqrtMagic = Math.Sqrt(magic);
            dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * Math.PI);
            dLon = (dLon * 180.0) / (a / sqrtMagic * Math.Cos(radLat) * Math.PI);
            double mgLat = lat + dLat;
            double mgLon = lon + dLon;           
            return new[] {mgLon, mgLat };		
		}

		/// <summary>
		/// Converts coordinates in projected meters to decimal degrees.
		/// </summary>
		/// <param name="p">Point in meters</param>
		/// <returns>Transformed point in decimal degrees</returns>
        protected override double[] MetersToRadians(double[] p)
		{
            double x = p[0] / 20037508.34 * 180;
            double y = p[1] / 20037508.34 * 180;
            y = 180.0 / Math.PI * (2 * Math.Atan(Math.Exp(y * Math.PI / 180.0)) - Math.PI / 2);
            return new double[] { x, y };       
		}
			
		/// <summary>
		/// Returns the inverse of this projection.
		/// </summary>
		/// <returns>IMathTransform that is the reverse of the current projection.</returns>
		public override IMathTransform Inverse()
		{
			if (_inverse==null)
                _inverse = new Wgs84toGcj02Math();
			return _inverse;
		}
    
    }

    public class Wgs84toGcj02Math : CoordinateSystems.Projections.MapProjection
    {
        public Wgs84toGcj02Math()
            : base(null)
        {
        }
   
        private double transformLat(double x, double y)
        {
            double ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * Math.Sqrt(Math.Abs(x));
            ret += (20.0 * Math.Sin(6.0 * x * Math.PI) + 20.0 * Math.Sin(2.0 * x * Math.PI)) * 2.0 / 3.0;
            ret += (20.0 * Math.Sin(y * Math.PI) + 40.0 * Math.Sin(y / 3.0 * Math.PI)) * 2.0 / 3.0;
            ret += (160.0 * Math.Sin(y / 12.0 * Math.PI) + 320 * Math.Sin(y * Math.PI / 30.0)) * 2.0 / 3.0;
            return ret;
        }

        private double transformLon(double x, double y)
        {
            double ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * Math.Sqrt(Math.Abs(x));
            ret += (20.0 * Math.Sin(6.0 * x * Math.PI) + 20.0 * Math.Sin(2.0 * x * Math.PI)) * 2.0 / 3.0;
            ret += (20.0 * Math.Sin(x * Math.PI) + 40.0 * Math.Sin(x / 3.0 * Math.PI)) * 2.0 / 3.0;
            ret += (150.0 * Math.Sin(x / 12.0 * Math.PI) + 300.0 * Math.Sin(x / 30.0 * Math.PI)) * 2.0 / 3.0;
            return ret;
        }
			
		/// <summary>
		/// Converts coordinates in decimal degrees to projected meters.
		/// </summary>
		/// <param name="lonlat">The point in decimal degrees.</param>
		/// <returns>Point in projected meters</returns>
        protected override double[] RadiansToMeters(double[] lonlat)
		{
            var lon = lonlat[0];		   
            var lat = lonlat[1];    

            double ee = 0.00669342162296594323;
            double a = 6378245.0;
            double dLat = transformLat(lon - 105.0, lat - 35.0);
            double dLon = transformLon(lon - 105.0, lat - 35.0);
            double radLat = lat / 180.0 * Math.PI;
            double magic = Math.Sin(radLat);
            magic = 1 - ee * magic * magic;
            double sqrtMagic = Math.Sqrt(magic);
            dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * Math.PI);
            dLon = (dLon * 180.0) / (a / sqrtMagic * Math.Cos(radLat) * Math.PI);
            double mgLat = lat + dLat;
            double mgLon = lon + dLon;           
            return new[] {mgLon, mgLat };		
		}

		/// <summary>
		/// Converts coordinates in projected meters to decimal degrees.
		/// </summary>
		/// <param name="p">Point in meters</param>
		/// <returns>Transformed point in decimal degrees</returns>
        protected override double[] MetersToRadians(double[] p)
		{
            double x = p[0] / 20037508.34 * 180;
            double y = p[1] / 20037508.34 * 180;
            y = 180.0 / Math.PI * (2 * Math.Atan(Math.Exp(y * Math.PI / 180.0)) - Math.PI / 2);
            return new double[] { x, y };       
		}
			
		/// <summary>
		/// Returns the inverse of this projection.
		/// </summary>
		/// <returns>IMathTransform that is the reverse of the current projection.</returns>
		public override IMathTransform Inverse()
		{
            if (_inverse == null)
                _inverse = new Gcj02toWgs84Math();
			return _inverse;
		}
    }

}
