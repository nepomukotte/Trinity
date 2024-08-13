import org.junit.* ;
import static org.junit.Assert.* ;
import CCfg.*;

public class JDocumentTest {
   @Test
   public void test_getConfig() {
      System.out.println("getConfig() Test") ;
      DocumentFacade myDoc = DocumentFacade.loadDoc("ConfigTest.xml");  
	  CConfigFacade conf = myDoc.getConfig();
	  assertEquals(conf.getName(""), "Test");
   }
   
	@Test
   public void test_saveToFile() {
      System.out.println("saveToFile() Test") ;
      DocumentFacade myDoc = DocumentFacade.loadDoc("ConfigTest.xml");  
	  CConfigFacade conf = myDoc.getConfig();
	  conf.setInt("a", "42");
	  myDoc.saveToFile("SavedConfigJava.xml");
	  myDoc = DocumentFacade.loadDoc("SavedConfigJava.xml");
	  conf = myDoc.getConfig();
	  assertEquals(conf.getInt("a"), 42); 
   }
	   
	@Test
   public void test_saveToUrl() {
      System.out.println("saveToUrl() Test") ;
      DocumentFacade myDoc = DocumentFacade.loadDoc("ConfigTest.xml");  
	  CConfigFacade conf = myDoc.getConfig();
	  conf.setInt("a", "42");
	  myDoc.saveToUrl("SavedConfigJava1.xml");
	  myDoc = DocumentFacade.loadDoc("SavedConfigJava1.xml");
	  conf = myDoc.getConfig();
	  assertEquals(conf.getInt("a"), 42); 
   }
	   
	@Test
   public void test_save() {
      System.out.println("save() Test") ;
      DocumentFacade mydoc = DocumentFacade.loadDoc("ConfigTest.xml");
      CConfigFacade conf = mydoc.getConfig();
      conf.removeNode("Numbers");
      mydoc.saveToFile("SavedConfigJava2.xml");
      mydoc = DocumentFacade.loadDoc("SavedConfigJava2.xml");
      conf = mydoc.getConfig();
      assertEquals(conf.testNode("Numbers"), false);
   }
}