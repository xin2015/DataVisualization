using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using System.Web.Mvc;

namespace DataVisualization.Web.Controllers
{
    public class OpenLayersController : Controller
    {
        // GET: OpenLayers
        public ActionResult Index()
        {
            return View();
        }

        public ActionResult Index2()
        {
            return View();
        }

        public ActionResult Index3()
        {
            return View();
        }

        public ActionResult Test()
        {
            return View();
        }

        #region Samples
        public ActionResult GeoJson()
        {
            return View();
        }

        public ActionResult CanvasGradientPattern()
        {
            return View();
        }
        #endregion
    }
}