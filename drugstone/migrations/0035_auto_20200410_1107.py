# Generated by Django 3.0.5 on 2020-04-10 11:07

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0034_auto_20200408_0952'),
    ]

    operations = [
        migrations.AlterField(
            model_name='drug',
            name='name',
            field=models.CharField(default='', max_length=256),
        ),
    ]