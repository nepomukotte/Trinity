import org.junit.* ;
import static org.junit.Assert.* ;
import CCfg.*;

public class JCConfigTest {
	 @Test
	   public void test_append() {
	      System.out.println("append() Test") ; 
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      conf.removeNode("");
	      assertTrue(conf.testNode("Numbers") == false) ;
	 } 
	 @Test
	   public void test_subConf() {
		 System.out.println("subConf() Test") ;       
		 DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
		 CConfigFacade conf = mydoc.getConfig();
		 CConfigFacade conf2 = conf.subConf("Numbers");
		 assertEquals(conf2.getName(""), "Numbers");
	 }
	 @Test
	   public void test_testNode() {
	      System.out.println("testNode() Test") ;
	      DocumentFacade  mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      assertEquals(conf.testNode("Numbers"), true);

	   }
	 @Test
	   public void test_createNode() {
	      System.out.println("createNode() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.createNode("NewNode");
	      assertEquals(conf.testNode("NewNode"), true);
	   }
	 @Test
	   public void test_removeNode() {
	      System.out.println("removeNode() Test") ;      
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.removeNode("Numbers");
	      assertEquals(conf.testNode("Numbers"), false);
	   }
	 @Test
	   public void test_enumNode() {
	      StringVector nodeList = new StringVector();
	      
	      System.out.println("enumNode() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.enumNode("", nodeList);
	      assertEquals(nodeList.get(0), "Numbers");
	   }
	 @Test
	   public void test_getCurrentPath() {
	      System.out.println("getCurrentPath() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      assertEquals(conf.getCurrentPath(), "Numbers");;
	   }
	 @Test
	   public void test_getInt() {
	      System.out.println("getInt() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      assertEquals(conf.getInt("i"), 3);
	   }
	 @Test
	   public void test_getStr() {
	      System.out.println("getStr() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setStr("str", "yop");
	      assertEquals(conf.getStr("str"), "yop");
	   } 
	 @Test
	   public void test_getBool() {
	      System.out.println("getBool() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setBool("bool", "true");
	      assertEquals(conf.getBool("bool"), true);
	   }
	 @Test
	   public void test_getReal() {
	      System.out.println("getReal() Test") ;
	      float temp = 4.2f;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setReal("real", "4.2");
	      assertEquals(conf.getReal("real"), temp, 0.001);
	   }
	 @Test
	   public void test_getName() {
	      System.out.println("getName() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      assertEquals(conf.getName(""), "Test");
	   }
	 @Test
	   public void test_getUnit() {
	      System.out.println("getUnit() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      assertEquals(conf.getUnit("memSize"), "bytes");
	   }
	 @Test
	   public void test_getIndex() {
	      System.out.println("getIndex() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      assertEquals(conf.getIndex("ind[0]"), "0");
	   }
	 @Test
	   public void test_getRange() {
	      System.out.println("getRange() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setRange("test", "[0, 100]");
	      assertEquals(conf.getRange("test"), "[0, 100]");
	   }
	 @Test
	   public void test_getPreComment() {
	      System.out.println("getPreComment() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setPreComment("test", "PreCom");
	      assertEquals(conf.getPreComment("test"), "PreCom");
	   }
	 @Test
	   public void test_getPostComment() {
	      System.out.println("getPostComment() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setPostComment("test", "PostCom");
	      assertEquals(conf.getPostComment("test"), "PostCom");
	   }
	 @Test
	   public void test_setCurrentPath() {
	      System.out.println("setCurrentPath() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setCurrentPath("Numbers");
	      assertEquals(conf.getCurrentPath(), "Numbers");
	   }
	 @Test
	   public void test_setInt() {
	      System.out.println("setInt() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      conf.setInt("j", "5");
	      assertEquals(conf.getInt("j"), 5);
	   }
	 @Test
	   public void test_setIntInt() {
	      System.out.println("setInt(int) Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      conf.setInt("j", 7);
	      assertEquals(conf.getAsStr("j"), "7");
	   }
	 @Test
	   public void test_setStr() {
	      System.out.println("setStr() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setStr("str", "yop");
	      assertEquals(conf.getStr("str"), "yop");
	   }
	 @Test
	   public void test_setBool() {
	      System.out.println("setBool() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setBool("bool", "true");
	      assertEquals(conf.getBool("bool"), true);
	   }
	 @Test
	   public void test_setBoolBoolean() {
	      System.out.println("setBool(boolean) Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setBool("bool", false);
	      assertEquals("false", conf.getAsStr("bool"));
	   }
	 @Test
	   public void test_setReal() {
	      System.out.println("setReal(String) Test") ;
	      float temp = 4.2f;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setReal("real", "4.2");
	      assertEquals(conf.getReal("real"), temp, 0.001f);
	   }
	 @Test
	   public void test_setRealFloat() {
	      System.out.println("setReal(float) Test") ;
	      float temp = 333.f;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setReal("real", temp);
	      assertEquals("333.0", conf.getAsStr("real"));
	   }
	 @Test
	   public void test_getAsInt() {
	      System.out.println("getAsInt() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      assertEquals(65536, conf.getAsInt("memSize"));
	   }
	 @Test
	   public void test_getHex() {
	      System.out.println("getHex() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      assertEquals(65536, conf.getHex("memSize"));
	   }
	 @Test
	   public void test_getHex32() {
	      System.out.println("getHex() 32-bit Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("hardwareDescription_fullCoBoStandAlone.xcfg");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Device[debug]");
	      assertEquals(2168717312L, conf.getHex("baseAddress"));
	   }
	 @Test
	   public void test_getAsStr() {
	      System.out.println("getAsStr() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("hardwareDescription_fullCoBoStandAlone.xcfg");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Device[debug]");
	      assertEquals("0x81440000", conf.getAsStr("baseAddress"));
	   }
	 @Test
	   public void test_setHexString() {
	      System.out.println("setHex(string) Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("hardwareDescription_fullCoBoStandAlone.xcfg");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Device[debug]");
	      conf.setHex("baseAddress", "0x100000");
	      assertEquals(1048576, conf.getHex("baseAddress"));
	   }
	 @Test
	   public void test_setHexLong() {
	      System.out.println("setHex(long) Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("hardwareDescription_fullCoBoStandAlone.xcfg");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Device[debug]");
	      conf.setHex("baseAddress", 0xF0F0F0F0);
	      assertEquals("0xF0F0F0F0", conf.getAsStr("baseAddress"));
	   }
	 @Test
	   public void test_setName() {
	      System.out.println("setName() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setName("", "Name");
	      assertEquals(conf.getName(""), "Name");
	   }
	 @Test
	   public void test_setUnit() {
	      System.out.println("setUnit() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.append("Numbers");
	      conf.setUnit("memSize", "meters");
	      assertEquals(conf.getUnit("memSize"), "meters");
	   }
	 @Test
	   public void test_setIndex() {
	      System.out.println("setIndex() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setIndex("ind[0]", "1");
	      assertEquals(conf.getIndex("ind[1]"), "1");
	   }
	 @Test
	   public void test_setRange() {
	      System.out.println("setRange() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setRange("test", "[0, 100]");
	      assertEquals(conf.getRange("test"), "[0, 100]");
	   }
	 @Test
	   public void test_setPreComment() {
	      System.out.println("setPreComment() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setPreComment("test", "PreCom");
	      assertEquals(conf.getPreComment("test"), "PreCom");
	   }
	 @Test
	   public void test_setPostComment() {
	      System.out.println("setPostComment() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      conf.setPostComment("test", "PostCom");
	      assertEquals(conf.getPostComment("test"), "PostCom");
	   }
	 @Test
	   public void test_getPathInt() {
	      System.out.println("getPathInt() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
	      CConfigFacade conf = mydoc.getConfig();
	      assertEquals(conf.getReal("Numbers.pos[x]"), 5.2, 1e-3);
	   }
	 @Test
	   public void test_getAgetShapingTime() {
	      System.out.println("getAgetShapingTime() Test") ;
	      DocumentFacade mydoc = DocumentFacade.loadDoc("configure-pedestals.xcfg");
	      CConfigFacade conf = mydoc.getConfig();
	      int shapingTime_ns = conf.getInt("Node[CoBo].Instance[0].AsAd[0].Aget[0].Global.Reg1.peackingTime");
	      System.out.println("Shaping time: " + shapingTime_ns + " ns");
	      boolean isSelected = conf.getBool("Node[CoBo].Instance[0].AsAd[0].Aget[0].channel[3].isSelectedforTestMode");
	      System.out.println("isSelected: " + isSelected);
	   }
}