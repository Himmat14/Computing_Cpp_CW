#include <iostream>
using std::string;

class AbstractEmployee {
    virtual void AskForPromotion()=0;
};

//blueprint for class
class Employee:AbstractEmployee {

public:
    //define attributes private access modifiers
    string Name;
    string Company;
    int Age;
    // attributes (private by default)
public:
    //setter and getter functions
    void setName(string name){
        Name = name;
    }
    string getName() {
        return Name;
    }

    void setAge(int age){
        if (age>=18)
            Age = age;
    }
    int getAge() {
        return Age;
    }

    void setCompany(string company){
        Company = company;
    }
    string getCompany() {
        return Company;
    }



    void IntroduceYourself() {
        std::cout << "Name - " << Name << std::endl;
        std::cout << "Company - " << Company << std::endl;
        std::cout << "Age - " << Age << std::endl;
    }
    //define Constructors
    /* constuctor must be public
     * no return type
     * same name as class*/
    Employee(string name, string company, int age) {
        Name = name;
        Company = company;
        Age = age;
    }

    void AskForPromotion() {
        if (Age>30) {
            std::cout << Name << " got promoted!!" << std::endl;
        }
        else {
            std::cout<<Name<< ",sorry NO promotion !"<<std::endl;
        }
    }

    virtual void Work() {
        std::cout<<Name<<" is checking email, task backlog, performing tasks..."<<std::endl;

    }

};

// developer clas is a child class of employee(parent)
class Developer: public Employee {
public:
    string FavProgrammingLanguage;

    Developer(string name, string company, int age, string favProgrammingLanguage)
            : Employee(name, company, age) {
        FavProgrammingLanguage = favProgrammingLanguage;
    }

    void FixBug() {
        std::cout << getName() << " fixed bug using " << FavProgrammingLanguage << std::endl;
    }


    void Work() {
        std::cout << Name << " is writting " << FavProgrammingLanguage << " code" << std::endl;

    };

};

class Teacher:public Employee {
public:
    string Subject;

    void PrepareLesson() {
        std::cout << Name << " is preparing " << Subject << " lesson" << std::endl;
    }

    //constructor
    Teacher(string name, string company, int age, string subject)
            : Employee(name, company, age) {
        Subject = subject;
    }

    void Work() {
        std::cout << Name << " is teaching " << Subject << std::endl;

    };

};

int main() {

    //define a emplioyee claass named emoployee1
    Employee employee1 = Employee("Himmy", "Imperial", 19);
    employee1.IntroduceYourself();


    Employee employee2 = Employee("Simmy", "Imperial", 69);
    employee2.IntroduceYourself();

    employee1.setAge(19);
    std::cout << employee1.getName() << " is " << employee1.getAge() << " years old." << std::endl;

    employee1.AskForPromotion();
    employee2.AskForPromotion();

    Developer d = Developer("himmy", "Imperial", 19, "C++");
    Teacher T = Teacher("simmy", "Imperial", 69, "French");

    T.PrepareLesson();
    T.AskForPromotion();
    d.FixBug();
    d.AskForPromotion();
    d.Work();
    T.Work();

    //pointers

    Employee* e1 = &d;
    Employee* e2 = &T;

    e1->Work();
    e2->Work();

};