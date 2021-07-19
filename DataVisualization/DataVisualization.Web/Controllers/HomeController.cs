using DataVisualization.Web.DataCenterService;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using System.Web.Mvc;

namespace DataVisualization.Web.Controllers
{
    public class HomeController : Controller
    {
        public ActionResult Index()
        {
            return View();
        }

        public ActionResult About()
        {
            ViewBag.Message = "Your application description page.";

            return View();
        }

        public ActionResult Contact()
        {
            ViewBag.Message = "Your contact page.";

            return View();
        }

        public ActionResult GetStationHourData()
        {
            using (DataCenterServiceClient client = new DataCenterServiceClient())
            {
                StationHourData[] stationHourData = client.GetStationHourDataListFromHistoryByTime("Suncere9999", "NULL4%tuNC3dO2fQ", DateTime.Today);
                return Json(stationHourData, JsonRequestBehavior.AllowGet);
            }
        }

        public ActionResult GetStationHourDataByTime(DateTime time)
        {
            using (DataCenterServiceClient client = new DataCenterServiceClient())
            {
                StationHourData[] stationHourData = client.GetStationHourDataListFromHistoryByTime("Suncere9999", "NULL4%tuNC3dO2fQ", time);
                return Json(stationHourData, JsonRequestBehavior.AllowGet);
            }
        }
    }
}